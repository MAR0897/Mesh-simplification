//=============================================================================
//
//  CLASS ModLindTurk - IMPLEMENTATION
//
//=============================================================================
#define OPENMESH_DECIMATER_MODLINDTURK_CC
//== INCLUDES =================================================================
#include <OpenMesh/Tools/Decimater/ModLindTurkT.hh>
//== NAMESPACE ===============================================================
namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Decimater { // BEGIN_NS_DECIMATER
//== IMPLEMENTATION ==========================================================

using Matrix3d = std::array<Vec3d, 3>;

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
initialize()
{
    //formalni inicializace vsech moznych veci, probiha pouze jednou, error zde nepocitej

    if (!LTprops.is_valid()) Base::mesh().add_property(LTprops);

    typename Mesh::HalfedgeIter he_it = Base::mesh().halfedges_begin(),
                                he_end = Base::mesh().halfedges_end();

    //Lock all boundary edges if option is set, lock all boundary
    // and "semi-boundary" edges, so that the mesh boundary stays the same
    if (lock_boundary_edges) {
        for (; he_it != he_end; ++he_it) {

            //no error was calculated yet
            Base::mesh().property(LTprops, *he_it).error_calculated = false;

            if (Base::mesh().is_boundary(*he_it)) {
                typename Mesh::VertexHandle vh1 = Base::mesh().to_vertex_handle(*he_it),
                                            vh2 = Base::mesh().from_vertex_handle(*he_it);
                typename Mesh::VertexOHalfedgeIter  voh_it1 = Base::mesh().voh_iter(vh1),
                                                    voh_it2 = Base::mesh().voh_iter(vh2);
                typename Mesh::VertexIHalfedgeIter  vih_it1 = Base::mesh().vih_iter(vh1),
                                                    vih_it2 = Base::mesh().vih_iter(vh2);
                for (; voh_it1.is_valid(); ++voh_it1) Base::mesh().property(LTprops, *voh_it1).is_locked = true;
                for (; voh_it2.is_valid(); ++voh_it2) Base::mesh().property(LTprops, *voh_it2).is_locked = true;
                for (; vih_it1.is_valid(); ++vih_it1) Base::mesh().property(LTprops, *vih_it1).is_locked = true;
                for (; vih_it2.is_valid(); ++vih_it2) Base::mesh().property(LTprops, *vih_it2).is_locked = true;
            }
        }
    }

    else for (; he_it != he_end; ++he_it) Base::mesh().property(LTprops, *he_it).error_calculated = false;

    std::cout<<"Options: "<<"\n\t"<<"Boundary locked? "<<std::boolalpha<<lock_boundary_edges<<"\n\t"<<
    "Lambda (error caculation weight): "<<"\t"<<lambda<<"\n\t"<<"Alpha parameter: "<<"\t"<<alpha<<std::endl;
}

template<class DecimaterType>
float
ModLindTurkT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci) 
{
    VertexHandle vh0 = _ci.v0;        //vertex to be potentially removed
    VertexHandle vh1 = _ci.v1;        //potentially remaining vertex
    HalfedgeHandle heh = _ci.v0v1;    //halfedge we are currently calculating error for

    if(!Base::mesh().property(LTprops, heh).is_locked and !Base::mesh().property(LTprops, _ci.v1v0).error_calculated){
        
        //set variables to zero
        Base::mesh().property(LTprops, heh).n = 0;                      
        Base::mesh().property(LTprops, heh).constraints.setZero();
        Base::mesh().property(LTprops, heh).b_side.setZero();           

        Eigen::Vector3d constraint; constraint.setZero();    //storage for new constraint
        double bside = 0.0;                                 //storage for new right side number

        //sets to avoid repeating calculations
        std::set<VertexHandle> vertex_handles;
        std::set<FaceHandle> face_handles;
        std::vector<Eigen::Vector3d> normals;
        std::vector<double> determinants;
           
        Eigen::Matrix3d Hv; Hv.setZero();           //Hessian for volume optimization
        Eigen::Matrix3d Hb; Hb.setZero();           //Hessian for boundary optimization
        Eigen::Matrix3d Hs; Hs.setZero();           //Hessian for triangle shape optimization
        Eigen::Vector3d cv; cv.setZero();           //vector for volume optimizaton
        Eigen::Vector3d cb; cb.setZero();           //vector for boundary optimization
        Eigen::Vector3d cs; cs.setZero();           //vector for triangle shape optimization
        double kv = 0.0;                            //constants in volume optimization
        double kb = 0.0;                            //constants in boundary optimization
        double ks = 0.0;                            //constants in triangle shape optimization
        Eigen::Matrix3d E1;  E1.setZero();          //e1 for every vertex
        Eigen::Matrix3d E2;  E2.setZero();          //e2 for every vertex
        Eigen::Matrix3d e1x; e1x.setZero();         //e1x matrix for boundary optimization
        Eigen::Vector3d e1; e1.setZero();                //summed E1        
        Eigen::Vector3d e2; e2.setZero();                //summed E2    
        Eigen::Vector3d e3; e3.setZero();                //cross product of e1 and e2
        size_t N = 0;                               //number of boundary edges if current edge is semiboundary (2 or 3)
        DefaultTraits::Point p;

        Eigen::Vector3d tri_shape; tri_shape.setZero();  //vector for storing vertex coords in triangle shape optimization    

    //---------------------------------------------------------------------------------------------------------------
    //Volume preservation (+volume optimization)    
        //plug the faces we need into a set (every element is unique)
        for (auto& vertex_handle : {vh0, vh1}) {
            typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(vertex_handle);
            for (; vf_it.is_valid(); ++vf_it) {
                face_handles.insert(*vf_it);
            } 
        }
        //calc constraint
        for (const auto& ff : face_handles) {

            //get the determinant of the face and compute the first bside
            Eigen::Matrix3d fv_coords;
            typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(ff);
            for (size_t i = 0; fv_it.is_valid(); ++fv_it, ++i) {
                DefaultTraits::Point p = Base::mesh().point(*fv_it);
                fv_coords.col(i) = Eigen::Vector3d(p[0], p[1], p[2]);
            }

            Eigen::Vector3d AB = fv_coords.col(1)-fv_coords.col(0);
            Eigen::Vector3d AC = fv_coords.col(2)-fv_coords.col(0);
            Eigen::Vector3d normal = AB.cross(AC);
            double determinant = fv_coords.col(0).dot(normal);

            constraint += normal;
            bside += determinant;
            
            //store values for possible calculation of remaining constraints after boundary preservation step
            normals.emplace_back(normal);
            determinants.emplace_back(determinant);
        }
        if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);

    //------------------------------------------------------------------------------------------------------------------------------
    //Boundary preservation
        if(Base::mesh().is_boundary(vh0) or Base::mesh().is_boundary(vh1)){
            //if edge is boundary, there will be 3 edges needed for constraints calculation (Figure 3), if not, there will be only 2
            if (Base::mesh().is_boundary(heh)) N = 3;
            else N = 2;
            std::vector<HalfedgeHandle> boundary_edges;

            if (N == 3) {
                HalfedgeHandle heh1 = Base::mesh().next_halfedge_handle(heh);
                HalfedgeHandle heh2 = Base::mesh().prev_halfedge_handle(heh);
                boundary_edges.emplace_back(heh1); 
                boundary_edges.emplace_back(heh2); 
                boundary_edges.emplace_back(heh);
            }
            else {
                //getting the two boundary edges
                typename Mesh::VertexOHalfedgeIter voh_it;
                typename Mesh::VertexIHalfedgeIter vih_it;
                //if v0 is boundary, check its outgoing halfedges, if not, check v0 outgoing halfedges
                if (Base::mesh().is_boundary(vh0)) {voh_it = Base::mesh().voh_iter(vh0); vih_it = Base::mesh().vih_iter(vh0);}
                else {voh_it = Base::mesh().voh_iter(vh1); vih_it = Base::mesh().vih_iter(vh1);}
                for (; voh_it.is_valid(); ++voh_it) if (Base::mesh().is_boundary(*voh_it)) boundary_edges.emplace_back(*voh_it);
                for (; vih_it.is_valid(); ++vih_it) if (Base::mesh().is_boundary(*vih_it)) boundary_edges.emplace_back(*vih_it);
            }

            //calculate E1 and E2 for every edge (that is for 2 or 3 edges)
            for (size_t i = 0; i<N; ++i) {
                HalfedgeHandle hh = boundary_edges[i];        
                VertexHandle vhto = Base::mesh().to_vertex_handle(hh);           
                VertexHandle vhfrom = Base::mesh().from_vertex_handle(hh);
                p = Base::mesh().point(vhto);   Eigen::Vector3d coords0  = Eigen::Vector3d(p[0], p[1], p[2]);
                p = Base::mesh().point(vhfrom); Eigen::Vector3d coords1  = Eigen::Vector3d(p[0], p[1], p[2]);
                E1.row(i) = coords1-coords0;
                E2.row(i) = coords1.cross(coords0);
                e1 += E1.row(i);
                e2 += E2.row(i);
            }
            e3 = e1.cross(e2);
            
            //equation 7
            constraint = e3*(e1.transpose()*e1);   bside = -(e3.transpose()*e3).value();
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);
            //equation 8
            constraint = e1.cross(e3);      bside = 0.0;
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);         
        }

    //-----------------------------------------------------------------------------------------------------------------------     
    //Volume optimization
        if(Base::mesh().property(LTprops, heh).n < 3) {
           
            //compute hessian, c and k using calculated values from volume preservation section
            size_t size = normals.size();
            for (size_t i = 0; i<size; ++i) {
                Hv += normals[i]*normals[i].transpose();
                cv -= determinants[i]*normals[i].transpose();
                kv += determinants[i]*determinants[i];
            }

            //rescale VertexOptimization variables to match the equation (9)
            Hv /= 18.0; cv /= 18.0; kv /= 18.0;

            calc_remaining_constraints(heh, Hv, cv);      
        }                  
    
    //----------------------------------------------------------------------------------------------------------------------
    //Boundary optimization
        if((Base::mesh().is_boundary(vh0) or Base::mesh().is_boundary(vh1)) and Base::mesh().property(LTprops, heh).n < 3) {

            //calculate hessian, c and k from values computed at the boundary preservation section
            for (size_t i = 0; i<N; ++i) {
                //create the (e x ) matrices
                e1x(0,1) = -E1(i,2); e1x(0,2) = E1(i,1); e1x(1,0) = E1(i,2);
                e1x(1,2) = -E1(i,0); e1x(2,0) = -E1(i,1); e1x(2,1) = E1(i,0);
                Hb += e1x*e1x.transpose();
                cb += (E1.row(i)).cross(E2.row(i));
                kb += (E2.row(i)*E2.row(i).transpose()).value();
            }

            //rescale BoundaryOptimization variables to match the equation (10)
            Hb *= 0.5; cb *= 0.5; kb *= 0.5;

            calc_remaining_constraints(heh, Hb, cb);
        }   

    //----------------------------------------------------------------------------------------------------------------------
    //Apply triangle shape opt. if necessary
        if(Base::mesh().property(LTprops, heh).n < 3){
            //insert needed vertices into a set
            for (auto& vv : {vh0, vh1}) {
                typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vv);
                for (; vv_it.is_valid(); ++vv_it) vertex_handles.insert(*vv_it);
            }
            //and erase those, which are not needed
            vertex_handles.erase(vh0);
            vertex_handles.erase(vh1);
            //calculate the Hessian and cs
            for (auto& vv : vertex_handles){
                p = Base::mesh().point(vv); tri_shape  = Eigen::Vector3d(p[0], p[1], p[2]);
                Hs(0,0)+=2; Hs(1,1)+=2; Hs(2,2)+=2;   //add identity matrix
                cs -= 2*tri_shape;
                ks += 2*(tri_shape.transpose()*tri_shape).value();
            }
            calc_remaining_constraints(heh, Hs, cs);
        }

    //----------------------------------------------------------------------------------------------------------------------
    //Calculate edge collapse error
        if(Base::mesh().property(LTprops, heh).n == 3){
            //get final vertex position (solve system of equations using inverse matrix)
            Eigen::Matrix3d A_inv = Base::mesh().property(LTprops, heh).constraints.inverse();
            Eigen::Vector3d bside = Base::mesh().property(LTprops, heh).b_side;
            Eigen::Vector3d V = A_inv*bside;
            //store the vertex position as Point
            for (int i = 0; i<3; ++i) Base::mesh().property(LTprops, heh).res_vertex_coords[i] = V[i];  

            //rescale VertexOptimization variables to match the equation (9)
            //Hv /= 18.0; cv /= 18.0; kv /= 18.0;
            //rescale BoundaryOptimization variables to match the equation (10)
            //Hb *= 0.5; cb *= 0.5; kb *= 0.5;
        
            //compute volume and boundary cost
            double fv = (0.5*(V.transpose()*(Hv*V)) + (cv.transpose()*V)).value() + 0.5*kv;  //volume objective function
            double fb = (0.5*(V.transpose()*(Hb*V)) + (cb.transpose()*V)).value() + 0.5*kb;  //area objective function
            //calculate final error
            p = Base::mesh().point(vh0); Eigen::Vector3d v0  = Eigen::Vector3d(p[0], p[1], p[2]);
            p = Base::mesh().point(vh1); Eigen::Vector3d v1  = Eigen::Vector3d(p[0], p[1], p[2]);
            double length = (v1-v0).norm();
            double err = lambda*fv +                    //volume opt
                        (1-lambda)*length*length*fb;    //boundary opt

            Base::mesh().property(LTprops, _ci.v0v1).error_calculated = true;

            return static_cast<float>(err); 
        } 
    }
    //if 2 or less contraints were found, we won't collapse this edge
    return FLT_MAX;
}

template<class DecimaterType>
bool
ModLindTurkT<DecimaterType>::
is_alpha_compatible(const HalfedgeHandle& heh, const Eigen::Vector3d& constraint)
{
    if(Base::mesh().property(LTprops, heh).n == 0){
        return  !(constraint(0) == 0.0 and constraint(1) == 0.0 and constraint(2) == 0.0);
    }
    else if(Base::mesh().property(LTprops, heh).n == 1){
        Eigen::Vector3d constraint0 = Base::mesh().property(LTprops, heh).constraints.row(0);
        return (std::pow(constraint0.transpose()*constraint, 2)
            < std::pow((constraint0.norm()*constraint.norm()),2)*COSALPHA2);
    }
    else if(Base::mesh().property(LTprops, heh).n == 2){
        Eigen::Vector3d crossp =  Base::mesh().property(LTprops, heh).constraints.row(0).cross(Base::mesh().property(LTprops, heh).constraints.row(1));
        return (std::pow((crossp.transpose()*constraint).value(), 2) 
            > std::pow((crossp.norm()*constraint.norm()),2)*SINALPHA2);
    }
    return false;
}

//Adds constraint to the system
template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
add_constraint(const HalfedgeHandle& heh, const Eigen::Vector3d& constraint, const double& right_side)
{
    Base::mesh().property(LTprops, heh).constraints.row(Base::mesh().property(LTprops, heh).n) = constraint;
    Base::mesh().property(LTprops, heh).b_side[Base::mesh().property(LTprops, heh).n] = right_side;
    Base::mesh().property(LTprops, heh).n++;
}

//Calculates remaining constraints if there are less than 3
template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
calc_remaining_constraints(const HalfedgeHandle& heh, const Eigen::Matrix3d& Hessian, const Eigen::Vector3d& c)
{
    size_t n = Base::mesh().property(LTprops, heh).n;
    size_t N = 3-n;
    //Create identity matrix of a size (3-n, 3)
        Eigen::MatrixXd I(N, 3);
        for (size_t i = 2-n, j = 2; i!=0; i--, j--) I(i,j) = 1;
    //Create orthogonal matrix Z
        Eigen::Matrix3d Z;
        Z = Base::mesh().property(LTprops, heh).constraints.transpose();
        if(n == 0) Z = Eigen::MatrixXd::Identity(3,3);             //If no constraints so far, create a matrix of standard base vectors
        else {
            if (n == 1) {Z(0,1) = Z(1,0); Z(1,1) = -Z(0,0);}    //Add first orthogonal vector
            Z.col(2) = Z.col(0).cross(Z.col(1));              //Add second orthogonal vector
        }
    //compute remaining constraints and b sides        
        auto temp = I*Z.inverse();
        auto constraints = temp*Hessian;
        auto bsides = -temp*c;
    //add constraints if possible
        for (size_t i = 0; i<=2-n; ++i)
            if(is_alpha_compatible(heh, constraints.row(i)))       
                add_constraint(heh, constraints.row(i), bsides(i));
}

// template<class DecimaterType>
// Matrix3d
// ModLindTurkT<DecimaterType>::inverse3x3(const Matrix3d& m){
//     double  a = m[0][0], b = m[1][0], c = m[2][0],
//             d = m[0][1], e = m[1][1], f = m[2][1],
//             g = m[0][2], h = m[1][2], i = m[2][2];
//     double  a00 = e*i-f*h,  a01 = -d*i+f*g, a02 = -g*e+d*h,
//             a10 = -b*i+c*h, a11 = a*i-c*g,  a12 = -a*h+b*g,
//             a20 = -e*c+f*b, a21 = -a*f+c*d, a22 = a*e-d*b;
//     Matrix3d adjugate;
//     adjugate[0] = Vec3d(a00, a01, a02);
//     adjugate[1] = Vec3d(a10, a11, a12);
//     adjugate[2] = Vec3d(a20, a21, a22);
//     double det = a*a00+b*a01+c*a02;
//     Matrix3d res;
//     for (int i = 0; i<3; ++i) res[i] = adjugate[i]*(1.0/det);
//     return res;
// }
//-----------------------------------------------------------------------------

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
