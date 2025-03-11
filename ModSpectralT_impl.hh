//=============================================================================
//
//  CLASS ModLindTurk - IMPLEMENTATION
//
//=============================================================================
#define OPENMESH_DECIMATER_MODSPECTRAL_CC
//== INCLUDES =================================================================
#include <OpenMesh/Tools/Decimater/ModSpectralT.hh>
#include <limits>
//== NAMESPACE ================================================================
namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Decimater { // BEGIN_NS_DECIMATER
//== IMPLEMENTATION ===========================================================

template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
initialize()
{
    //-------------------------------------------------------------------------
    // OBECNE VECI
    //-------------------------------------------------------------------------

    //add props
    if (!SPprops.is_valid()) Base::mesh().add_property(SPprops);
    if (!idx.is_valid()) Base::mesh().add_property(idx);
    if (!local_idx.is_valid()) Base::mesh().add_property(local_idx);
    if (!inner.is_valid()) Base::mesh().add_property(inner);
    if (!area.is_valid()) Base::mesh().add_property(area);
    if (!cotangents.is_valid()) Base::mesh().add_property(cotangents);

    //Lock all boundary edges if option is set, lock all boundary
    // and "semi-boundary" edges, so that the mesh boundary stays the same
    typename Mesh::HalfedgeIter he_it = Base::mesh().halfedges_begin(),
                                he_end = Base::mesh().halfedges_end();
    if (lock_boundary_edges) {
        for (; he_it != he_end; ++he_it) {
            if (Base::mesh().is_boundary(*he_it)) {
                typename Mesh::VertexHandle                         //define halfedge's vertices
                    vh1 = Base::mesh().to_vertex_handle(*he_it),
                    vh2 = Base::mesh().from_vertex_handle(*he_it);
                typename Mesh::VertexOHalfedgeIter                  //get the vertices' outgoing halfedges
                    voh_it1 = Base::mesh().voh_iter(vh1),
                    voh_it2 = Base::mesh().voh_iter(vh2);
                typename Mesh::VertexIHalfedgeIter                  //get the vertices' ingoing halfedges
                    vih_it1 = Base::mesh().vih_iter(vh1),
                    vih_it2 = Base::mesh().vih_iter(vh2);
                for (; voh_it1.is_valid(); ++voh_it1)               //set the lock parameter to true for all
                    Base::mesh().property(SPprops, *voh_it1).is_locked = true;
                for (; voh_it2.is_valid(); ++voh_it2) 
                    Base::mesh().property(SPprops, *voh_it2).is_locked = true;
                for (; vih_it1.is_valid(); ++vih_it1) 
                    Base::mesh().property(SPprops, *vih_it1).is_locked = true;
                for (; vih_it2.is_valid(); ++vih_it2) 
                    Base::mesh().property(SPprops, *vih_it2).is_locked = true;
            }
        }
    }

    he_it = Base::mesh().halfedges_begin();
    for (; he_it!=he_end; ++he_it) Base::mesh().property(SPprops, *he_it).error_calculated = false;

    //-------------------------------------------------------------------------
    // VECI PRO SPEKTRALNI SIMPLIFIKACI
    //-------------------------------------------------------------------------

    size_t n_vertices = Base::mesh().n_vertices();          //number of vertices in the mesh
    norms.setZero(n_vertices);                              //previous error
    Eigen::VectorXd M;  M.setZero(n_vertices);              //declare mass matrix
    Eigen::SparseMatrix<double> L(n_vertices, n_vertices);  //declare laplacian matrix
    L.setZero();


    //MASS matrix calculation (loop through all vertices and sum up their neighboring
    // faces areas) and INDEXING the vertices for easier vector and matrix access
    typename Mesh::VertexIter v_it = Base::mesh().vertices_begin(),
                              v_end = Base::mesh().vertices_end();
    for (size_t i = 0; v_it!=v_end; ++v_it, ++i) {                          //iterate over the vertices
        Base::mesh().property(idx, *v_it) = i;                              //INDEX the vertices
        Base::mesh().property(inner, *v_it) = false;                        //no vertex is currently in 1R circle of a edge
        typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(*v_it);
        for (; vf_it.is_valid(); ++vf_it) {                                 //get the neighboring faces
            double obsah = Base::mesh().calc_face_area(*vf_it);             //calculate the face area
            M[i] += obsah;                                                  //add it to the vertex index in the mass matrix
            Base::mesh().property(area, *vf_it) = obsah;                    //store it as a property to reduce the amount of calculations 
        }
        M[i] /= 3.0;                                                          //divide by 3 (see mass matrix definition)
    }

    //LAPLACIAN matrix of the mesh
    v_it = Base::mesh().vertices_begin();
    for (; v_it!=v_end; ++v_it) {                                           //iterate over the mesh vertices
        size_t i = Base::mesh().property(idx, *v_it);                       //get the vertex index
        double sum = 0.0;                                                   //sum for laplacian matrix diagonal elements
        typename Mesh::VertexEdgeIter ve_it = Base::mesh().ve_iter(*v_it);                                                        
        for (; ve_it.is_valid(); ++ve_it) {                                 //iterate over 1 ring edges (cant use voh_it, cause it wont work on boundary edges)
            HalfedgeHandle he1 = Base::mesh().halfedge_handle(*ve_it, 0);                               //get first halfedge
            HalfedgeHandle he2 = Base::mesh().halfedge_handle(*ve_it, 1);                               //get second halfedge
            VertexHandle vh2;                                                                           //get 1-ring vertex
            if (Base::mesh().from_vertex_handle(he1) != *v_it) vh2 = Base::mesh().from_vertex_handle(he1); 
            else vh2 = Base::mesh().from_vertex_handle(he2);
            size_t j = Base::mesh().property(idx, vh2);                                                 //get its index
            if (!Base::mesh().is_boundary(he1)) {                                                   //calc cotangent alpha
                double cotg = /*0.5**/calc_cotangent(he1, *v_it, vh2);                                      //get the cotangent                
                sum += cotg;                                                                                //add it to the sum variable
                L.coeffRef(i,j) -= cotg;                                                                    //add it to its place in the laplacian matrix
            } 
            if (!Base::mesh().is_boundary(he2)) {                                                       //cotangent beta
                double cotg = /*0.5**/calc_cotangent(he2, *v_it, vh2);                                          //get the cotangent
                sum += cotg;                                                                                //add it to the sum variable
                L.coeffRef(i,j) -= cotg;                                                                    //add it to its place in the laplacian matrix
            }
        }
        L.coeffRef(i,i) = sum;                                                                     //add the minus sum to the diagonal
    }

    //get first eigenvec_n smallest eigenvectors of the laplacian matrix (= F matrix (signals))
	Eigen::VectorXd N = M.cwiseSqrt().cwiseInverse();
	Eigen::SparseMatrix<double> H = N.asDiagonal() * L * N.asDiagonal();
	size_t ncv = std::min(3 * eigenvec_n, static_cast<int>(L.rows()));
	Spectra::SparseSymShiftSolve<double> op(H);
    //Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> solver(op, eigenvec_n, ncv, -1e-6);
    Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, decltype(op)> solver(&op, eigenvec_n, ncv, -1e-6);
	solver.init();
	//solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestAlge);
	solver.compute(1000, 1e-10, Spectra::SMALLEST_ALGE);

	F = N.asDiagonal() * solver.eigenvectors();
	F.colwise().normalize();                    //signals
    //compute Z matrix
	Z = M.cwiseInverse().asDiagonal() * L * F;  //signals_L
	F.transposeInPlace();
	Z.transposeInPlace();

    Eigen::VectorXd LL;  LL.setZero(n_vertices); 
	Z.row(0) = LL;

    // std::cerr << "Mass: \n" <<  M << std::endl;
    // std::cerr << "Laplacian: \n" << L << std::endl;
    // std::cerr << "Signals: \n" <<  F << std::endl;
    // std::cerr << "Signals_L: \n" << Z << std::endl;

    // typename Mesh::FaceIter f_it = Base::mesh().faces_begin(),
    //                         f_end = Base::mesh().faces_end();
    // std::cerr << "AREAS:" << std::endl;
    // for (; f_it!=f_end; ++f_it) std::cerr << Base::mesh().property(area, *f_it) << "\n";
	// std::cerr << "COTANGENTS" << std::endl;
    // f_it = Base::mesh().faces_begin();
	// for (; f_it!=f_end; ++f_it) {
    //     int i = 0;
    //     for (auto& el : Base::mesh().property(cotangents, *f_it)) {
    //         if (i++ < 2) std::cerr << el.second << "\t";
    //         else std::cerr << el.second;
    //     }
    //     std::cerr << "\n";
    // }
	// std::cerr << std::endl;


    std::cout<<"Successfully initialized"<<std::endl;

    int NN = 0;
    typename Mesh::EdgeIter e_it = Base::mesh().edges_begin(),
                             e_end = Base::mesh().edges_end();
    for (; e_it!=e_end; ++e_it, ++NN); 

    std::cout <<"Number of edges = " << NN << std::endl; 
}
    
int GGGcount = 0;
//=============================================================================

template<class DecimaterType>
float
ModSpectralT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci) 
{
    HalfedgeHandle heh = _ci.v0v1;    //halfedge we are currently calculating error for

    if(!Base::mesh().property(SPprops, _ci.v0v1).is_locked /*and !Base::mesh().property(SPprops, _ci.v1v0).error_calculated*/){
        
        //--------------------------------------------------------------------------------------
        // Clear previous local indexing
        //--------------------------------------------------------------------------------------
        typename Mesh::VertexIter v_it = Base::mesh().vertices_begin(),
                                  v_end = Base::mesh().vertices_end();
        for(; v_it!=v_end; ++v_it) Base::mesh().property(local_idx, *v_it) = -1;

        //--------------------------------------------------------------------------------------
        // Get needed sets and index all LOCAL (2-ring) vertices for easies matrices access
        //--------------------------------------------------------------------------------------
        std::set<VertexHandle> verts1R; //1-ring vertices (around the edge)     //innerV in the original
        std::set<VertexHandle> verts2R; //2-ring of vertices (around the edge)  //V in the original
        std::set<FaceHandle> faces2R;   //2-ring of faces
        for (VertexHandle vertex : {_ci.v0, _ci.v1}) {
            typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vertex);
            for (; vv_it.is_valid(); ++vv_it) {                                            
                verts1R.insert(*vv_it);
                typename Mesh::VertexVertexIter vvv_it = Base::mesh().vv_iter(*vv_it);
                for (; vvv_it.is_valid(); ++vvv_it) verts2R.insert(*vvv_it);
                typename Mesh::VertexFaceIter vvf_it = Base::mesh().vf_iter(*vv_it);
                for (; vvf_it.is_valid(); ++vvf_it) faces2R.insert(*vvf_it);
            }
        }

        // std::cerr << "verts1R/innerV = " << verts1R.size() << "\t"
        //           << "verts2R/V = " << verts2R.size() << "\t"
        //           << "faces2R/F = " << faces2R.size() << "\t" <<std::endl;
       
        //---------------------------------------------------------------------
        // INDEX all LOCAL (2ring) vertices
        //---------------------------------------------------------------------
        int i = 0;
        for (auto& v : verts1R) Base::mesh().property(local_idx, v) = i++;
        for (auto& v : verts2R) if (Base::mesh().property(local_idx, v) == -1) Base::mesh().property(local_idx, v) = i++; 
        int v1R_size = verts1R.size(); 
        //add up the previous local cost (ERROR) from norms
        double prev_cost_local = 0;
        //std::cerr << "normy: ";
        for(auto& v : verts1R) {
            prev_cost_local += norms[Base::mesh().property(idx, v)];
            //std::cerr << norms[Base::mesh().property(idx, v)] << "; ";
        }
        
        //std::cerr << std::endl;

        //---------------------------------------------------------------------
        //setting PF matrix (precompute signals restriction)
        //---------------------------------------------------------------------
        Eigen::MatrixXd PF(F.rows(), verts2R.size());                        //a copy of signals matrix (F) only for local vertices
        Eigen::SparseMatrix<double> L(verts2R.size(), verts2R.size());    //modified laplacian matrix (L s vlnovkou)
        Eigen::MatrixXd LPF(F.rows(), verts2R.size());                       //PF*L
        for(auto& vertex : verts2R) {                                        //get the PF matrix
            size_t index = Base::mesh().property(idx, vertex);
            int local_index = Base::mesh().property(local_idx, vertex);
            if(vertex == _ci.v0) PF.col(local_index).setConstant(std::numeric_limits<double>::quiet_NaN());
            else if(vertex != _ci.v1) PF.col(local_index) = F.col(index);
        }

        //std::cerr << "Signal restriction matrix: \n" << PF << std::endl;


        //---------------------------------------------------------------------
        // Function that evaluates collapse error for different value of where
        // on the edge the resulting vertex will be (from 0 to 1). In this algorithm
        // the value are always 0.0, 0.5 and 1.0, which are then used to calculate
        // minimum of a polymonial function f : <0,1> -> R+_0. For the minimum,
        // this function will then also calculate the error (if not between some tolerances)
        //---------------------------------------------------------------------
        auto eval = [&](double alpha, std::vector<std::pair<VertexHandle, double>>& diff) -> std::pair<Vec3d, double>
        {
            Vec3d p0 = vector_cast<Vec3d>(Base::mesh().point(_ci.v0));      //remove vertex point (_ci.v0)
            Vec3d p1 = vector_cast<Vec3d>(Base::mesh().point(_ci.v1));      //keep vertex point (_ci.v1)
            size_t idxR = Base::mesh().property(idx, _ci.v0);               //remove vertex index (_ci.v0)
            size_t idxK = Base::mesh().property(idx, _ci.v1);               //keep vertex index (_ci.v1)
            int local_idxR = Base::mesh().property(local_idx, _ci.v0);      //remove vertex local index (_ci.v0)
            int local_idxK = Base::mesh().property(local_idx, _ci.v1);      //keep vertex local index (_ci.v1)
            const Vec3d pos = p1*(1-alpha) + p0*alpha;                      //final position on the edge (this function goes for alpha = 0, 0.5 and 1 to minimize cost polynomial)

            //restrict signals
            PF.col(local_idxK) = (F.col(idxK) * (1 - alpha) + F.col(idxR) * alpha);

            //std::cerr << "PF col: " << PF.col(local_idxK) << "\t";

            //set up thing for modified mass and laplacian matrix
            Eigen::VectorXd M; M.setZero(verts2R.size());
            std::vector<Eigen::Triplet<double>> coeffs; coeffs.clear();

            int IIII = 0;

            //go through all 2ring faces and exclude the 2 triangles that will get removed by the edge collapse
            for(auto& f : faces2R) if ((Base::mesh().face_handle(heh) != f and 
                                       Base::mesh().face_handle(Base::mesh().opposite_halfedge_handle(heh)) != f ) and
                                       !Base::mesh().status(f).deleted() ) { //mozna neni treba
                IIII++;
                //check if the face has _ci.v0 or _ci.v1 and get the vertices for calculations
                bool has_central_vertex = false;
                std::vector<VertexHandle> faceV;    //vertices of one face
                for (typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(f); fv_it.is_valid(); ++fv_it) {
                    if (*fv_it == _ci.v0 or *fv_it == _ci.v1) has_central_vertex = true;  //if we find that a face has one of the collapsing vertices, we need to calculate it separately
                    else faceV.emplace_back(*fv_it);    //and we will also collect the two other vertices
                }

                // temp vars for storing cotangents and area
                std::unordered_map<VertexHandle, double> c = Base::mesh().property(cotangents, f);
                double a = Base::mesh().property(area, f);

                // Function for adding values to the modified Laplacian
                auto add_coeffs = [&local_idxR, &coeffs, &v1R_size](int row, int col, double val){
                    if(col == local_idxR or col >= v1R_size) return;    // if the column vertex is the remove vertex or is not in 1-ring dont do anything
                    coeffs.emplace_back(row, col, val);
                };

                //recalculate the contangets and area if the faces are going to be deformed
                if(has_central_vertex){
                    VertexHandle v1 = faceV[0];
                    VertexHandle v2 = faceV[1];
                    Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(v1));
                    Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(v2));
                    Vec3d vec3 = pos;
                    Vec3d AB = vec1-vec3;
                    Vec3d AC = vec2-vec3;
                    a = 0.5*(AB.cross(AC)).norm();                              //calculate new (modified face) area
                    c.clear();
                    c.insert(std::pair<VertexHandle, double>(_ci.v1, calc_cotangent_from_points(vec1, vec2, vec3))); //lets assume the KEEP (_ci.v1) vertex has the final pos
                    c.insert(std::pair<VertexHandle, double>(v1, calc_cotangent_from_points(vec3, vec2, vec1)));
                    c.insert(std::pair<VertexHandle, double>(v2, calc_cotangent_from_points(vec1, vec3, vec2)));
                    for (auto& v : {v1, v2, _ci.v1}) M[Base::mesh().property(local_idx, v)] += a/3;    //add mass to modified mass matrix
                    std::vector<size_t> indexes;
                    std::vector<double> cotgs;
                    for (auto it = c.begin(); it!=c.end(); ++it){
                        indexes.emplace_back(Base::mesh().property(local_idx, it->first));
                        cotgs.emplace_back(it->second);
                    }   
                    // Laplacian
                    for (size_t k = 0; k<3; ++k) {                                                  //add cotangets to modified laplacian matrix  
                        add_coeffs(indexes[k], indexes[k], /*-0.5**/(cotgs[(k+1)%3]+cotgs[(k+2)%3]));
                        add_coeffs(indexes[k], indexes[(k+1)%3], -/*0.5**/cotgs[(k+2)%3]);
                        add_coeffs(indexes[k], indexes[(k+2)%3], -/*0.5**/cotgs[(k+1)%3]);
                    }
                }
                
                // if nothing needs to be recalculed, simply add stuff to the M and L
                else {
                    for (auto& v : faceV) M[Base::mesh().property(local_idx, v)] += a/3;            //add mass to modified mass matrix
                    std::vector<size_t> indexes;
                    std::vector<double> cotgs;
                    for (auto it = c.begin(); it!=c.end(); ++it){
                        indexes.emplace_back(Base::mesh().property(local_idx, it->first));
                        cotgs.emplace_back(it->second);
                    }
                    for (size_t k = 0; k<3; ++k) {                                                  //add cotangets to modified laplacian matrix 
                        add_coeffs(indexes[k], indexes[k], /*0.5**/(cotgs[(k+1)%3]+cotgs[(k+2)%3]));
                        add_coeffs(indexes[k], indexes[(k+1)%3], -/*0.5**/cotgs[(k+2)%3]);
                        add_coeffs(indexes[k], indexes[(k+2)%3], -/*0.5**/cotgs[(k+1)%3]);
                    }
                }
            }

            
            //set up the modified laplacian matrix and compute error
            L.setFromTriplets(coeffs.begin(), coeffs.end());

            //std::cerr << L << std::endl;
            LPF = PF * L;
            //compute local cost (sum of E_v) (minus the previous local cost)
            double cost_local = 0.0;
            for(auto& v : verts1R) if(v != _ci.v0) {
                size_t l_idx = Base::mesh().property(local_idx, v);                                 //get local index
                size_t g_idx = Base::mesh().property(idx, v);                                       //get global index of same vertex
                double Mv = M[l_idx];                                                               //take out mass for that vertex
                double E_v = std::numeric_limits<double>::quiet_NaN();                              //initialize error with maximal value
                if(v == _ci.v1){       
                    const Eigen::VectorXd PZv = Z.col(idxK) * (1 - alpha) + Z.col(idxR) * alpha;
                    E_v = Mv * (PZv - (1.0 / Mv) * LPF.col(l_idx)).squaredNorm();                   //compute squared norm for _ci.v1
                }
                else E_v = Mv * (Z.col(g_idx) - (1.0 / Mv) * LPF.col(l_idx)).squaredNorm();         //compute squared norm
                cost_local += E_v;                                                                  //sum up E_v
                diff.emplace_back(v, E_v);                                                          //store for previous cost for next collapse edge 
            }
            const double cost = cost_local - 0.75*prev_cost_local;
            // std::cerr << "["<<IIII<<"]"<< "/" << Base::mesh().property(idx, _ci.v1) << "-" << Base::mesh().property(idx, _ci.v0)  << "/" <<"For "<<alpha<<": " << cost << std::endl;
            // std::cerr << "COST LOCAL bez norem: " << cost_local << std::endl;
            return {pos, cost};
        };


        //---------------------------------------------------------------------
        // calculate the error, minimize it, update halfedge properties and return
        //---------------------------------------------------------------------
        std::vector<std::pair<VertexHandle, double>> diff[3];
        std::pair<Vec3d, double> cost[3] = { eval(0.0, diff[0]), eval(0.5, diff[1]), eval(1.0, diff[2]) };  
        
        //minimize polynomial
        double y1 = cost[0].second;
        double y2 = cost[1].second;
        double y3 = cost[2].second;
        double minimum = (-(-3*y1+4*y2-y3)) / (2*(2*y1-4*y2+2*y3));            // -b/2a

        //dont calculate if not needed, else recalculate for the computed minimum
        std::pair<Vec3d, double> final_edge_cost;
        std::vector<std::pair<VertexHandle, double>> cd;
        if (minimum < 0.005) { minimum = 0.0; final_edge_cost = cost[0]; cd = diff[0]; }                       //if minimum is < 0.005, use the calculated value for 0.0
        else if (std::abs(minimum-0.5) < 0.005) { minimum = 0.5; final_edge_cost = cost[1]; cd = diff[1]; }    //if minimum is around 0.5, use the calculated value for 0.5
        else if (minimum > 0.995) { minimum = 1.0; final_edge_cost = cost[2]; cd = diff[2]; }                  //if minimum is > 0.995, use the calculated value for 1.0
        else final_edge_cost = eval(minimum, cd);                                                              //else calculate it again for the minimum

        //mozna budes muset to minimum overit, zdali je skutecne minimum, bo nekdy to muze vyjit i s vetsim errorem (parabola je naopak)
        double min_y = std::min({y1, y2, y3});
        if (final_edge_cost.second > min_y) {

            if (y1 <= y2 and y1 <= y3) { minimum = 0.0; final_edge_cost = cost[0]; cd = diff[0]; }
            else if (y3 <= y2) { minimum = 1.0; final_edge_cost = cost[2]; cd = diff[2]; }
            else { minimum = 0.5; final_edge_cost = cost[1]; cd = diff[1]; }
        }

        //update mesh properties
        Base::mesh().property(SPprops, heh).alpha = minimum;                            
        Base::mesh().property(SPprops, heh).cost_diff = cd;                             //store cost_diff for norms, which is then used for previous local cost in next collapse
        Base::mesh().property(SPprops, heh).res_vertex_coords = final_edge_cost.first;  //resulting collapse vertex coordinates
        Base::mesh().property(SPprops, heh).error_calculated = true;
            Base::mesh().property(SPprops, heh).error = final_edge_cost.second;
            Base::mesh().property(SPprops, heh).prev_local_cost = prev_cost_local;

        // double el = (Base::mesh().point(_ci.v0)-Base::mesh().point(_ci.v1)).norm();
        // std::cerr << "/"<< Base::mesh().property(idx, _ci.v0) << "-" <<Base::mesh().property(idx, _ci.v1) 
        // << "/"<< "Edge cost: " << final_edge_cost.second << "\t" << "Edge length: "<< el << "\t" << 
        // "Prev cost: " << prev_cost_local << std::endl;

        return static_cast<float>(final_edge_cost.second);   
    }
    
    return FLT_MAX;
}

//=============================================================================

template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
preprocess_collapse(const CollapseInfo& _ci)
{
    //move remaining vertex to ideal calculated position
    DefaultTraits::Point ideal_vertex;
    for (int i = 0; i<3; ++i) ideal_vertex[i] = 
        Base::mesh().property(SPprops, _ci.v0v1).res_vertex_coords[i];
    Base::mesh().set_point(_ci.v1, ideal_vertex);
}


template<class DecimaterType>
void
ModSpectralT<DecimaterType>::
postprocess_collapse(const CollapseInfo& _ci)
{


	// std::cerr << "{"<< ++GGGcount << "}"<< "/"<< Base::mesh().property(idx, _ci.v0)
    // << "-"<< Base::mesh().property(idx, _ci.v1)<< "/"<<"edge zkolabovan------------------------------------------------" << std::endl;
    // std::cerr << "Error kolapsu = " << Base::mesh().property(SPprops, _ci.v0v1).error << "\t" << "PrevLocalCost = " << Base::mesh().property(SPprops, _ci.v0v1).prev_local_cost << std::endl;
    // std::cerr << "----------------------------------------------" << std::endl;
	

    // typename Mesh::VertexIter v_it = Base::mesh().vertices_begin(),
    //                           v_end = Base::mesh().vertices_end();
    // std::cerr << "Vertex indexes: " << std::endl; 
    // for (; v_it != v_end; ++v_it) if (!Base::mesh().status(*v_it).deleted()) std::cerr << Base::mesh().property(idx, *v_it) << "\t";
    // std::cerr << std::endl;


    //mozna zkus to udelat pro opacny halfedge, ten prece zustava
    double alpha = Base::mesh().property(SPprops, _ci.v0v1).alpha;
    std::vector<std::pair<VertexHandle, double>> cost_diff = 
        Base::mesh().property(SPprops, _ci.v0v1).cost_diff;
    // Update signals
    size_t idxR = Base::mesh().property(idx, _ci.v0);
    size_t idxK = Base::mesh().property(idx, _ci.v1);
	F.col(idxK) = F.col(idxK) * (1 - alpha) + F.col(idxR) * alpha;
	Z.col(idxK) = Z.col(idxK) * (1 - alpha) + Z.col(idxR) * alpha;
	F.col(idxR).setConstant(std::numeric_limits<double>::quiet_NaN());
	Z.col(idxR).setConstant(std::numeric_limits<double>::quiet_NaN());

    // Update costs
	for(const std::pair<VertexHandle, double>& d : cost_diff) {
        norms[Base::mesh().property(idx, d.first)] = d.second;
        //std::cerr << "{idx" << Base::mesh().property(idx, d.first) << "}"<< "COST DIFF UPDATE " << d.second << "\n";
    }

    //std::cerr << "vektor NOREM: \n" << norms << std::endl;

    //get faces that need to recalculate
    std::vector<FaceHandle> faces;
    typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(_ci.v1);
    for (; vf_it.is_valid(); ++vf_it) faces.emplace_back(*vf_it);
    std::vector<VertexHandle> v;
    //recalculate each face
    for (auto& f : faces) {
        v.clear();
        typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(f);
        for (; fv_it.is_valid(); ++fv_it) v.emplace_back(*fv_it);
        Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(v[0]));
        Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(v[1]));
        Vec3d vec3 = vector_cast<Vec3d>(Base::mesh().point(v[2]));
        Vec3d AB = vec1-vec3;
        Vec3d AC = vec2-vec3;
        double a = 0.5*(AB.cross(AC)).norm();
        std::unordered_map<VertexHandle, double> cotgs;
        cotgs.insert(std::pair<VertexHandle, double>(v[2], calc_cotangent_from_points(vec1, vec2, vec3)));
        cotgs.insert(std::pair<VertexHandle, double>(v[0], calc_cotangent_from_points(vec3, vec2, vec1)));
        cotgs.insert(std::pair<VertexHandle, double>(v[1], calc_cotangent_from_points(vec1, vec3, vec2)));
        Base::mesh().property(area, f) = a;
        Base::mesh().property(cotangents, f) = cotgs;
    }


    // typename Mesh::FaceIter f_it = Base::mesh().faces_begin(),
    //                         f_end = Base::mesh().faces_end();
    // std::cerr << "Updated areas after collapse" << std::endl;
    // for (; f_it!=f_end; ++f_it) std::cerr << Base::mesh().property(area, *f_it) << "\n";
	// std::cerr << "Updated cotangets after collapse" << std::endl;
    // f_it = Base::mesh().faces_begin();
	// for (; f_it!=f_end; ++f_it) {
    //     int i = 0;
    //     for (auto& el : Base::mesh().property(cotangents, *f_it)) {
    //         if (i++ < 2) std::cerr << el.second << "\t";
    //         else std::cerr << el.second;
    //     }
    //     std::cerr << "\n";
    // }
	// std::cerr << std::endl;

}

//=============================================================================

template<typename T>
double 
ModSpectralT<T>::calc_cotangent(const HalfedgeHandle& he, const VertexHandle& vh1, const VertexHandle& vh2){
    FaceHandle fh = Base::mesh().face_handle(he);
    typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(fh);
    VertexHandle vh3;
    for (; fv_it.is_valid(); ++fv_it) if (*fv_it!=vh1 and *fv_it!=vh2) vh3 = *fv_it; //get 3rd vertex  
    Vec3d vec1 = vector_cast<Vec3d>(Base::mesh().point(vh1));
    Vec3d vec2 = vector_cast<Vec3d>(Base::mesh().point(vh2));
    Vec3d vec3 = vector_cast<Vec3d>(Base::mesh().point(vh3));
    double cotg = calc_cotangent_from_points(vec1, vec2, vec3);
    Base::mesh().property(cotangents, fh).insert(std::pair<VertexHandle, double>(vh3, cotg));
    return cotg;
}

//=============================================================================

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


//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
