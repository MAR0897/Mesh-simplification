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

    //std::cout<<"F: \n"<<F<<std::endl;

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
        std::set<FaceHandle> faces1R;       //1-ring of faces for recalculation of mass and laplacian
        for (VertexHandle vertex : {vh0, vh1}) {
            for (typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vertex); vv_it.is_valid(); ++vv_it) {
                vertices1R.insert(*vv_it);
                for (typename Mesh::VertexVertexIter vvv_it = Base::mesh().vv_iter(*vv_it); vvv_it.is_valid(); ++vvv_it) vertices2R.insert(*vvv_it);
                for (typename Mesh::VertexFaceIter vvf_it = Base::mesh().vf_iter(*vv_it); vvf_it.is_valid(); ++vvf_it) faces2R.insert(*vvf_it);
            }
            for (typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(vertex); vf_it.is_valid(); ++vf_it) faces1R.insert(*vf_it);
        }
        std::cout<<"segfault0"<<std::endl;
        // faces to then recalculate mass and cotangents for if the collapse is done
        for (auto it = faces1R.begin(); it != faces1R.end(); ) {
            if (Base::mesh().face_handle(heh) == *it || Base::mesh().face_handle(Base::mesh().opposite_halfedge_handle(heh)) == *it) {
                it = faces1R.erase(it);  // erase returns the next valid iterator
            } else {
                ++it;  // only increment if not erasing
            }
        }
        size_t i = 0;
        for (auto& v : vertices2R) Base::mesh().property(local_idx, v) = i++;  //index all local (2ring) vertices
        std::cout<<"segault0.1"<<std::endl;
        // Previous local cost
        double prev_cost_local = 0;
        for(auto& v : vertices1R) prev_cost_local += norms[Base::mesh().property(idx, v)];
            
        //setting PF matrix (Precompute signals restriction)
        Eigen::MatrixXd PF(F.rows(), vertices2R.size());
        Eigen::SparseMatrix<double> L(vertices2R.size(), vertices2R.size());
        Eigen::MatrixXd LPF(F.rows(), vertices2R.size());
        for(auto& vertex : vertices2R) {
            size_t index = Base::mesh().property(idx, vertex);
            size_t local_index = Base::mesh().property(local_idx, vertex);
            if(vertex == vh0) PF.col(local_index).setConstant(std::numeric_limits<double>::quiet_NaN());
            else if(vertex != vh1) PF.col(local_index) = F.col(index);
        }

        auto eval = [&](double alpha, std::vector<std::pair<VertexHandle, double>>& diff) -> std::pair<Vec3d, double>
        {
            Vec3d p0 = vector_cast<Vec3d>(Base::mesh().point(vh0)); //REMOVE
            Vec3d p1 = vector_cast<Vec3d>(Base::mesh().point(vh1)); //KEEP
            size_t idxR = Base::mesh().property(idx, vh0);          //remove vertex index
            size_t idxK = Base::mesh().property(idx, vh1);          //keep vertex index
            size_t local_idxR = Base::mesh().property(local_idx, vh0);          //remove vertex local index
            size_t local_idxK = Base::mesh().property(local_idx, vh1);          //keep vertex local index
            const Vec3d pos = p1*(1-alpha) + p0*alpha;              //final position on the edge (this function goes for alpha = 0, 0.5 and 1 to minimize cost polynomial)

            // Restrict signals
            PF.col(local_idxK) = (F.col(idxK) * (1 - alpha) + F.col(idxR) * alpha);

            Eigen::VectorXd M = Eigen::VectorXd::Zero(vertices2R.size());
            std::vector<Eigen::Triplet<double>> coeffs;
            coeffs.clear();

            std::cout<<"segfault 1"<<std::endl;
            //Go through all 2ring faces and exclude the 2 triangles that will get removed by the edge collapse
            for(auto& f : faces2R) if (Base::mesh().face_handle(heh) != f and Base::mesh().face_handle(Base::mesh().opposite_halfedge_handle(heh)) != f) {

                //check if the face vh0 or vh1 and get the vertices for calculations
                bool has_central_vertex = false;
                std::vector<VertexHandle> faceV;
                for (typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(f); fv_it.is_valid(); ++fv_it) {
                    if (*fv_it == vh0 or *fv_it == vh1) has_central_vertex = true;
                    else faceV.emplace_back(*fv_it);
                }
                std::cout<<"segfault 2"<<std::endl;
                //recalculate the contangets and area if the faces are going to be deformed
                std::unordered_map<VertexHandle, double> c = Base::mesh().property(cotangents, f);
                double a = Base::mesh().property(area, f);
                if(has_central_vertex){
                    VertexHandle V1 = faceV[0];
                    VertexHandle V2 = faceV[1];
                    Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(V1));
                    Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(V2));
                    Vec3d vec3 = pos;
                    Vec3d AB = vec1-vec3;
                    Vec3d AC = vec2-vec3;
                    a = 0.5*(AB.cross(AC)).norm();
                    c.clear();
                    c.insert(std::pair<VertexHandle, double>(vh1, calc_cotangent_from_points(vec1, vec2, vec3))); //lets assume the KEEP vertex has the final pos
                    c.insert(std::pair<VertexHandle, double>(V1, calc_cotangent_from_points(vec3, vec2, vec1)));
                    c.insert(std::pair<VertexHandle, double>(V2, calc_cotangent_from_points(vec1, vec3, vec2)));
                    std::cout<<"segfault 4"<<std::endl;
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
                        coeffs.emplace_back(indexes[k], indexes[k], -(cotgs[(k+1)%3]+cotgs[(k+2)%3]));//???
                        coeffs.emplace_back(indexes[k], indexes[(k+1)%3], 0.5*cotgs[(k+2)%3]);//???
                        coeffs.emplace_back(indexes[k], indexes[(k+2)%3], 0.5*cotgs[(k+1)%3]);//???
                    }
                }
                std::cout<<"segfault 3"<<std::endl;
                //std::cout<<"eval func: faces loop: M and L recalculated"<<std::endl;
            }
            std::cout<<"segfault 5"<<std::endl;
            if(iiiii > 14875) {
                std::cout<<vertices2R.size()<<std::endl;
                for (auto e : coeffs) std::cout<<e.row()<<", "<<e.col()<<", "<<e.value()<<std::endl;
            }
            L.setFromTriplets(coeffs.begin(), coeffs.end());
            std::cout<<"segfault 5.1"<<std::endl;
            // Cost (local) (sum of E_v) (minus the previous local cost)
            double cost_local = 0.0;
            LPF = PF * L;
            for(auto& v : vertices1R) if(v != vh0) {
                size_t idx = Base::mesh().property(local_idx, v);
                double Mv = M[idx];
                double E_v = std::numeric_limits<double>::quiet_NaN();
                if(v == vh1){
                    const Eigen::VectorXd PZv = Z.col(idxK) * (1 - alpha) + Z.col(idxR) * alpha;
                    E_v = Mv * (PZv - (1.0 / Mv) * LPF.col(idx)).squaredNorm();
                }
                else E_v = Mv * (Z.col(idx) - (1.0 / Mv) * LPF.col(idx)).squaredNorm();
                cost_local += E_v;
                diff.emplace_back(v, E_v);
            }
            std::cout<<"segfault 6"<<std::endl;
            const double cost = cost_local - prev_cost_local;
            return {pos, cost};
        };


        std::vector<std::pair<VertexHandle, double>> diff[4];
        std::pair<Vec3d, double> cost[3] = { eval(0.0, diff[0]), eval(0.5, diff[1]), eval(1.0, diff[2]) };

        //minimize polynom
        double y1 = cost[0].second;
        double y2 = cost[1].second;
        double y3 = cost[2].second;
        double minimum = -(-y3+4*y2-3*y1) / 2*(2*y3-4*y2+2*y1);            // -b/2a
        if (minimum<0.0) minimum = 0.0;
        else if (minimum>1.0) minimum = 1.0;

        std::pair<Vec3d, double> final_edge_cost = eval(minimum, diff[3]);       //TODO: set tolerances to lower the number of calculations (if the minimum is 0.0, we already calculated this value)

        Base::mesh().property(SPprops, heh).alpha = minimum;
        Base::mesh().property(SPprops, heh).cost_diff = diff[3];
        Base::mesh().property(SPprops, heh).res_vertex_coords = final_edge_cost.first;
        Base::mesh().property(SPprops, heh).recalc_faces = faces1R;
        return static_cast<float>(final_edge_cost.second);   
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
    //std::cout<<norms<<std::endl;

    //TODO: recalculate the cotangents
    std::set<FaceHandle> faces = Base::mesh().property(SPprops, _ci.v0v1).recalc_faces;
    for (auto& f : faces) {
        typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(f);
        VertexHandle v1 = *fv_it;
        VertexHandle v2 = *(++fv_it);
        VertexHandle v3 = *(++fv_it);
        Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(v1));
        Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(v2));
        Vec3d vec3 = vector_cast<Vec3d>(Base::mesh().point(v3));
        Vec3d AB = vec1-vec3;
        Vec3d AC = vec2-vec3;
        double a = 0.5*(AB.cross(AC)).norm();
        std::unordered_map<VertexHandle, double> cotgs;
        cotgs.insert(std::pair<VertexHandle, double>(v3, calc_cotangent_from_points(vec1, vec2, vec3))); //lets assume the KEEP vertex has the final pos
        cotgs.insert(std::pair<VertexHandle, double>(v1, calc_cotangent_from_points(vec3, vec2, vec1)));
        cotgs.insert(std::pair<VertexHandle, double>(v1, calc_cotangent_from_points(vec1, vec3, vec2)));

        Base::mesh().property(area, f) = a;
        Base::mesh().property(cotangents, f) = cotgs;
    }
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
