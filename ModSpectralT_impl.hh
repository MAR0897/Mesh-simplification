//=============================================================================
//
//  CLASS ModLindTurk - IMPLEMENTATION
//
//=============================================================================
#define OPENMESH_DECIMATER_MODSPECTRAL_CC
//== INCLUDES =================================================================
#include <OpenMesh/Tools/Decimater/ModSpectralT.hh>
#include <limits>
//== NAMESPACE ===============================================================
namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Decimater { // BEGIN_NS_DECIMATER
//== IMPLEMENTATION ==========================================================


template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
initialize()
{
    
    //formalni inicializace vsech moznych veci, probiha pouze jednou, error zde nepocitej
    //-------obecne veci
    if (!SPprops.is_valid()) Base::mesh().add_property(SPprops);
    if (!idx.is_valid()) Base::mesh().add_property(idx);
    if (!local_idx.is_valid()) Base::mesh().add_property(local_idx);
    if (!area.is_valid()) Base::mesh().add_property(area);
    if (!cotangents.is_valid()) Base::mesh().add_property(cotangents);

    typename Mesh::HalfedgeIter he_it = Base::mesh().halfedges_begin(),
                                he_end = Base::mesh().halfedges_end();
    //Lock all boundary edges if option is set, lock all boundary
    // and "semi-boundary" edges, so that the mesh boundary stays the same
    if (lock_boundary_edges) {
        for (; he_it != he_end; ++he_it) {
            if (Base::mesh().is_boundary(*he_it)) {
                typename Mesh::VertexHandle vh1 = Base::mesh().to_vertex_handle(*he_it),
                                            vh2 = Base::mesh().from_vertex_handle(*he_it);
                typename Mesh::VertexOHalfedgeIter  voh_it1 = Base::mesh().voh_iter(vh1),
                                                    voh_it2 = Base::mesh().voh_iter(vh2);
                typename Mesh::VertexIHalfedgeIter  vih_it1 = Base::mesh().vih_iter(vh1),
                                                    vih_it2 = Base::mesh().vih_iter(vh2);
                for (; voh_it1.is_valid(); ++voh_it1) Base::mesh().property(SPprops, *voh_it1).is_locked = true;
                for (; voh_it2.is_valid(); ++voh_it2) Base::mesh().property(SPprops, *voh_it2).is_locked = true;
                for (; vih_it1.is_valid(); ++vih_it1) Base::mesh().property(SPprops, *vih_it1).is_locked = true;
                for (; vih_it2.is_valid(); ++vih_it2) Base::mesh().property(SPprops, *vih_it2).is_locked = true;
            }
        }
    }

    //----------veci pro spektralni spimlifikaci
    size_t n_vertices = Base::mesh().n_vertices();
    norms.setZero(n_vertices);
	projection.resize(n_vertices, n_vertices);
	projection.setIdentity();

    //mass matrix calculation
    Eigen::VectorXd M = Eigen::VectorXd::Zero(n_vertices);
    typename Mesh::VertexIter v_it = Base::mesh().vertices_begin();
    typename Mesh::VertexIter v_end = Base::mesh().vertices_end();
    for (size_t i = 0; v_it!=v_end; ++v_it, ++i) {
        Base::mesh().property(idx, *v_it) = i;  //INDEX the vertices
        typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(*v_it);
        for (; vf_it.is_valid(); ++vf_it) {
            double obsah = Base::mesh().calc_face_area(*vf_it);
            M[i] += obsah;
            Base::mesh().property(area, *vf_it) = obsah;
        }
        M[i] /= 3;
    }

    //laplacian operator of the mesh
    L.resize(n_vertices, n_vertices);
    v_it = Base::mesh().vertices_begin();
    v_end = Base::mesh().vertices_end();
    for (; v_it!=v_end; ++v_it) {  //iterate over the whole mesh
        size_t i = Base::mesh().property(idx, *v_it); 
        typename Mesh::VertexOHalfedgeIter voh_it = Base::mesh().voh_iter(*v_it);
        double sum = 0.0;
        for (; voh_it.is_valid(); ++voh_it){ //iterate over the outgoing halfedges
            VertexHandle vh2 = Base::mesh().to_vertex_handle(*voh_it); //1-ring vertex
            size_t j = Base::mesh().property(idx, vh2); 
            HalfedgeHandle he2 = Base::mesh().opposite_halfedge_handle(*voh_it); //opposite halfedge to get the two bordering faces
            if (!Base::mesh().is_boundary(*voh_it)) {   //cotangent alpha
                double cotg = calc_cotangent(*voh_it, *v_it, vh2);
                sum += cotg;
                L.coeffRef(i,j) += 0.5*cotg; 
            } 
            if (!Base::mesh().is_boundary(he2)) {   //cotangent beta
                double cotg = calc_cotangent(he2, *v_it, vh2);
                sum += cotg;
                L.coeffRef(i,j) += 0.5*cotg; 
            }
        }
        L.coeffRef(i,i) = -sum;
    }
	
	const Eigen::VectorXd N = M.cwiseSqrt().cwiseInverse();
	const Eigen::SparseMatrix<double> H = N.asDiagonal() * L * N.asDiagonal();
	const unsigned int ncv = std::min(3 * eigenvec_n, static_cast<int>(L.rows()));
	Spectra::SparseSymShiftSolve<double> op(H);
    //u kralika nejde, LU faktorizace nebo neco selze
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> solver(op, eigenvec_n, ncv, -1e-6);
	solver.init();
	solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestAlge);
	F = N.asDiagonal() * solver.eigenvectors();
	F.colwise().normalize();
	Z = M.cwiseInverse().asDiagonal() * L * F;
	F.transposeInPlace();
	Z.transposeInPlace();

    std::cout<<"Successfully initialized"<<std::endl;
}
//=======================================================================================================================

template<class DecimaterType>
float
ModSpectralT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci) 
{
    VertexHandle vh0 = _ci.v0;        //vertex to be potentially removed
    VertexHandle vh1 = _ci.v1;        //potentially remaining vertex
    HalfedgeHandle heh = _ci.v0v1;    //halfedge we are currently calculating error for

    if(!Base::mesh().property(SPprops, heh).is_locked){

        // Get 2-rings and 1-ring (inner vertices??) (and get other iterable stuff)
        std::set<VertexHandle> vertices1R;  //1-ring vertices (around the edge)     //innerV in the original
        std::set<VertexHandle> vertices2R;  //2-ring of vertices (around the edge)  //V in the original
        std::set<FaceHandle> faces2R;       //2-ring of faces (correct??)
        for (VertexHandle vertex : {vh0, vh1}) {
            typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vh0);
            for (; vv_it.is_valid(); ++vv_it) {
                vertices1R.insert(*vv_it);
                typename Mesh::VertexVertexIter vvv_it = Base::mesh().vv_iter(*vv_it);
                for (; vvv_it.is_valid(); ++vvv_it) vertices2R.insert(*vvv_it);
                typename Mesh::VertexFaceIter vvf_it = Base::mesh().vf_iter(*vv_it);
                for (; vvf_it.is_valid(); ++vvf_it) faces2R.insert(*vvf_it);
            }
        }
        std::set<VertexHandle> vertices1Rout = vertices1R; //outer 1 ring of vertices (without vh0 and vh1 = to be removed and kept vertices)
        vertices1Rout.erase(vh0);
        vertices1Rout.erase(vh1);
        std::set<VertexHandle> vertices2Rout = vertices2R;  //outer 2 ring of vertices (without whole 1 ring)
        for (auto v : vertices1R) vertices2Rout.erase(v);

        std::cout<<"Two rings: gottem"<<std::endl;

        // Previous local cost
        double prev_cost_local = 0;
        for(auto& v : vertices1R) prev_cost_local += norms[Base::mesh().property(idx, v)];
        
        std::cout<<"Prev local costmaxxed"<<std::endl;

        //setting PF matrix (Precompute signals restriction)
        Eigen::MatrixXd PF(F.rows(), vertices2R.size());
        Eigen::SparseMatrix<double> L(vertices2R.size(), vertices2R.size());
        Eigen::MatrixXd LPF(F.rows(), vertices2R.size());
        for(VertexHandle vertex : vertices2R) {
            size_t index = Base::mesh().property(idx, vertex);
            if(vertex == vh0) PF.col(index).setConstant(std::numeric_limits<double>::quiet_NaN());
            else if(vertex != vh1) PF.col(index) = F.col(index);
        }

        std::cout<<"PF letter sent"<<std::endl;


        auto eval = [&](double alpha, std::vector<std::pair<VertexHandle, double>>& diff) -> std::pair<Vec3d, double>
        {
            Vec3d p0 = vector_cast<Vec3d>(Base::mesh().point(vh0)); //REMOVE
            Vec3d p1 = vector_cast<Vec3d>(Base::mesh().point(vh1)); //KEEP
            size_t idxR = Base::mesh().property(idx, vh0);          //remove vertex index
            size_t idxK = Base::mesh().property(idx, vh1);          //keep vertex index
            const Vec3d pos = p1*(1-alpha) + p0*alpha;              //final position on the edge (this function goes for alpha = 0, 0.5 and 1 to minimize cost polynomial)
            
            //if we have one of vh0 and vh1 we select the resulting position instead (makes sense)
            //auto new_pos = [&](uint32_t i) { return i == to_keep || i == to_remove ? pos : object->vertices[i]; };
            //auto new_index = [&](uint32_t i) { return i == to_remove ? to_keep : i; };

            // Restrict signals
            PF.col(idxK) = (F.col(idxK) * (1 - alpha) + F.col(idxR) * alpha);

            Eigen::VectorXd M = Eigen::VectorXd::Zero(vertices2R.size());
            size_t i = 0;
            for (auto& v : vertices2R) if (v != vh0) Base::mesh().property(local_idx, v) = i++;  //index all remaining vertices

            std::vector<Eigen::Triplet<double>> coeffs;
            coeffs.clear();
            for(auto& f : faces2R){
                // Exclude the 2 triangles that will get removed by the edge collapse
                if (Base::mesh().face_handle(heh) == f or Base::mesh().face_handle(Base::mesh().opposite_halfedge_handle(heh)) == f) continue;
                
                //check if the face vh0 or vh1 and get the vertices for calculations
                bool has_central_vertex = false;
                std::vector<VertexHandle> faceV;
                typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(f);
                for (; fv_it.is_valid(); ++fv_it) {
                    if (*fv_it == vh0 or *fv_it == vh1) has_central_vertex = true;
                    else faceV.emplace_back(*fv_it);
                }

                //recalculate the contangets and area if the faces are going to be deformed
                std::unordered_map<VertexHandle, double> c = Base::mesh().property(cotangents, f);
                double a = Base::mesh().property(area, f);
                if(has_central_vertex){
                    VertexHandle V1 = faceV[0];
                    VertexHandle V2 = faceV[1];
                    typename Mesh::Point p1 = Base::mesh().point(V1);
                    typename Mesh::Point p2 = Base::mesh().point(V2);
                    Vec3d vec1 = vector_cast<Vec3d>(p1);
                    Vec3d vec2 = vector_cast<Vec3d>(p2);
                    Vec3d vec3 = pos;
                    Vec3d AB = vec1-vec3;
                    Vec3d AC = vec2-vec3;
                    a = 0.5*(AB.cross(AC)).norm();
                    c.clear();
                    c.insert(std::pair<VertexHandle, double>(vh1, calc_cotangent_from_points(vec1, vec2, vec3))); //lets assume the KEEP vertex has the final pos
                    c.insert(std::pair<VertexHandle, double>(V1, calc_cotangent_from_points(vec3, vec2, vec1)));
                    c.insert(std::pair<VertexHandle, double>(V2, calc_cotangent_from_points(vec1, vec3, vec2)));

                    for (auto& v : {V1, V2, vh1}) M[Base::mesh().property(local_idx, v)] += a/3;
                    std::vector<size_t> indexes;
                    std::vector<double> cotgs;
                    for (auto it = c.begin(); it!=c.end(); ++it){
                        indexes.emplace_back(Base::mesh().property(local_idx, it->first));
                        cotgs.emplace_back(it->second);
                    }
                    for (size_t k = 0; k<3; ++k) {
                        coeffs.emplace_back(indexes[k], indexes[k], -cotgs[(k+1)%3]-cotgs[(k+2)%3]);//???
                        coeffs.emplace_back(indexes[k], indexes[(k+1)%3], -cotgs[(k+2)%3]);//???
                        coeffs.emplace_back(indexes[k], indexes[(k+2)%3], -cotgs[(k+1)%3]);//???
                    }
                }
                
                // if nothing needs to be recalculed, simply add stuff to the M and L
                else {
                    for (auto& v : faceV) M[Base::mesh().property(local_idx, v)] += a/3;
                    std::vector<size_t> indexes;
                    std::vector<double> cotgs;
                    for (auto it = c.begin(); it!=c.end(); ++it){
                        indexes.emplace_back(Base::mesh().property(local_idx, it->first));
                        cotgs.emplace_back(it->second);
                    }
                    for (size_t k = 0; k<3; ++k) {
                        coeffs.emplace_back(indexes[k], indexes[k], -cotgs[(k+1)%3]-cotgs[(k+2)%3]);//???
                        coeffs.emplace_back(indexes[k], indexes[(k+1)%3], -cotgs[(k+2)%3]);//???
                        coeffs.emplace_back(indexes[k], indexes[(k+2)%3], -cotgs[(k+1)%3]);//???
                    }
                }
            }
            L.setFromTriplets(coeffs.begin(), coeffs.end());

            // Cost (local) (sum of E_v) (minus the previous local cost)
            double cost_local = 0.0;
            LPF = PF * L;
            for(auto& v : vertices1R){
                if(v == vh0) continue;
                size_t idx = Base::mesh().property(local_idx, v);
                const double Mv = M[idx];
                double E_v = std::numeric_limits<double>::quiet_NaN();
                if(v == vh1){
                    const Eigen::VectorXd PZv = Z.col(idxK) * (1 - alpha) + Z.col(idxR) * alpha;
                    E_v = Mv * (PZv - (1.0 / Mv) * LPF.col(idx)).squaredNorm();
                }
                else E_v = Mv * (Z.col(idx) - (1.0 / Mv) * LPF.col(idx)).squaredNorm();
                cost_local += E_v;
                diff.emplace_back(v, E_v);
            }
            const double cost = cost_local - prev_cost_local;
            return {pos, cost};
        };


        std::vector<std::pair<VertexHandle, double>> diff[4];
        std::pair<Vec3d, double> cost[3] = { eval(0.0, diff[0]), eval(0.5, diff[1]), eval(1.0, diff[2]) };

        std::cout<<"Baseness evaluated at least once"<<std::endl;

        //minimize polynom
        double y1 = cost[0].second;
        double y2 = cost[1].second;
        double y3 = cost[2].second;
        double minimum = -(-y3+4*y2-3*y1) / 2*(2*y3-4*y2+2*y1);            // -b/2a
        if (minimum<0.0) minimum = 0.0;
        else if (minimum>1.0) minimum = 1.0;
        
        std::cout<<"Polynom destroyed"<<std::endl;

        std::pair<Vec3d, double> final_edge_cost = eval(minimum, diff[3]);       //TODO: set tolerances to lower the number of calculations (is the minimum is 0.0, we already calculated this value)
        
        Base::mesh().property(SPprops, heh).alpha = minimum;
        Base::mesh().property(SPprops, heh).cost_diff = diff[3];
        Base::mesh().property(SPprops, heh).res_vertex_coords = final_edge_cost.first;
        return static_cast<float>(final_edge_cost.second); 

        std::cout<<"Costcore finished"<<std::endl;
    }

    return FLT_MAX;
}

//=======================================================================================================================

template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
preprocess_collapse(const CollapseInfo& _ci)
{
    //move remaining vertex to ideal calculated position
    
    DefaultTraits::Point ideal_vertex;
    for (int i = 0; i<3; ++i) ideal_vertex[i] = Base::mesh().property(SPprops, _ci.v0v1).res_vertex_coords[i];
    Base::mesh().set_point(_ci.v1, ideal_vertex);
}


template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
postprocess_collapse(const CollapseInfo& _ci)
{

    //mozna zkus to udelat pro opacny halfedge, ten prece zustava
    double alpha = Base::mesh().property(SPprops, _ci.v0v1).alpha;
    std::vector<std::pair<VertexHandle, double>> cost_diff = Base::mesh().property(SPprops, _ci.v0v1).cost_diff;
    // Update signals
    size_t idxR = Base::mesh().property(idx, _ci.v0);
    size_t idxK = Base::mesh().property(idx, _ci.v1);
	F.col(idxK) = F.col(idxK) * (1 - alpha) + F.col(idxR) * alpha;
	Z.col(idxK) = Z.col(idxK) * (1 - alpha) + Z.col(idxR) * alpha;
	F.col(idxR).setConstant(std::numeric_limits<double>::quiet_NaN());
	Z.col(idxR).setConstant(std::numeric_limits<double>::quiet_NaN());

	// Update projection
	projection.row(idxK) = projection.row(idxK) * (1 - alpha) + projection.row(idxR) * alpha;
	//projection.prune([&](int row, int /*col*/, double /*val*/){ return row != (int)to_remove; });

    // Update costs
	for(const std::pair<VertexHandle, double>& d : cost_diff) norms[Base::mesh().property(idx, d.first)] = d.second;
	total_cost = norms.sum();
}

template<typename T>
double 
ModSpectralT<T>::calc_cotangent(const HalfedgeHandle& he, const VertexHandle& vh1, const VertexHandle& vh2){
    FaceHandle fh1 = Base::mesh().face_handle(he);
    typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(fh1);
    VertexHandle vh3;
    for (; fv_it.is_valid(); ++fv_it) if (*fv_it!=vh1 and *fv_it!=vh2) vh3 = *fv_it; //get 3rd vertex  
    Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(vh1));
    Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(vh2));
    Vec3d vec3 = vector_cast<Vec3d>(Base::mesh().point(vh3));
    double cotg = calc_cotangent_from_points(vec1, vec2, vec3);
    Base::mesh().property(cotangents, fh1).insert(std::pair<VertexHandle, double>(vh3, cotg));
    return cotg;
}

template<typename T>
double
ModSpectralT<T>::calc_cotangent_from_points(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3){
    Vec3d vec1 = p1-p3;
    Vec3d vec2 = p2-p3;
    double dot = (vec1.dot(vec2)) / (vec1.norm()*vec2.norm());
    double arccos = acos(dot); 
    double cotg = cos(arccos)/sin(arccos);
    return cotg;
}

//-----------------------------------------------------------------------------

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
