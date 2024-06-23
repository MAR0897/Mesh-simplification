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

    std::cout<<"Options: "<<lock_boundary_edges<<"\t"<<lambda<<"\t"<<alpha<<std::endl;
}

template<class DecimaterType>
float
ModLindTurkT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci) 
{
    VertexHandle vh0 = _ci.v0;        //vertex to be potentially removed
    VertexHandle vh1 = _ci.v1;        //potentially remaining vertex
    HalfedgeHandle heh = _ci.v0v1;    //halfedge we are currently calculating error for
    
    //HalfedgeHandle rev_heh = _ci.v1v0;   //reverse halfedgehandle
    /*if (A) Base::mesh().property(LTprops, heh).ncalc=0;
    Base::mesh().property(LTprops, heh).ncalc++;*/

    if(!Base::mesh().property(LTprops, heh).is_locked){
        
        //set variables to zero
        Vec3d zeros(0.0);   //Helper zeros-filler vector
        Base::mesh().property(LTprops, heh).n = 0;                      
        Base::mesh().property(LTprops, heh).constraints.fill(zeros);
        Base::mesh().property(LTprops, heh).b_side = zeros;           
        Base::mesh().property(LTprops, heh).res_vertex_coords = zeros;

        Vec3d constraint(0.0);                      //storage for new constraint
        double bside = 0.0;                         //storage for new right side number

        //sets to avoid repeating calculations
        std::set<VertexHandle> vertex_handles;
        std::set<FaceHandle> face_handles;
           
        Matrix3d Hv;            //Hessian for volume optimization
        Matrix3d Hb;            //Hessian for boundary optimization
        Matrix3d Hs;            //Hessian for triangle shape optimization
        Hv.fill(zeros);
        Hb.fill(zeros);
        Hs.fill(zeros);
        Vec3d cv(0.0);          //vector for volume optimizaton
        Vec3d cb(0.0);          //vector for boundary optimization
        Vec3d cs(0.0);          //vector for triangle shape optimization
        double kv = 0.0;        //constants in volume optimization
        double kb = 0.0;        //constants in boundary optimization
        double ks = 0.0;        //constants in triangle shape optimization

        Matrix3d E1;               //e1 for every vertex
        Matrix3d E2;              //e2 for every vertex
        Matrix3d e1x;               //e1x matrix for boundary optimization
        E1.fill(zeros);
        E2.fill(zeros);
        e1x.fill(zeros);              
        Vec3d e1(0.0);              //summed E1
        Vec3d e2(0.0);               //summed E2
        Vec3d e3(0.0);               //cross product of e1 and e2
        
        Vec3d tri_shape(0.0);       //vector for storing vertex coords in triangle shape optimization

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
        for (const auto& face_handle : face_handles) {
            if(Base::mesh().is_valid_handle(face_handle)){

                //get the determinant of the face and compute the first bside
                Matrix3d fv_coords;
                typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(face_handle);
                for (size_t i = 0; fv_it.is_valid(); ++fv_it, ++i) 
                    fv_coords[i] = vector_cast<Vec3d>(Base::mesh().point(*fv_it));
                    
                Vec3d AB = fv_coords[1]-fv_coords[0];
                Vec3d AC = fv_coords[2]-fv_coords[0];
                Vec3d normal = AB.cross(AC);
                constraint += normal;
                double determinant = fv_coords[0].dot(normal);
                bside += determinant;
                
                //Volume optimization section
                for (int i = 0; i<3; ++i) {
                    for (int j = 0; j<3; ++j) {
                        Hv[i][j] += normal[i]*normal[j];
                    }
                }
                cv -= determinant*normal;
                kv += determinant*determinant;
            }
            
        }
        if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);

    //------------------------------------------------------------------------------------------------------------------------------
    //Boundary preservation
        if(Base::mesh().is_boundary(vh0) or Base::mesh().is_boundary(vh1)){
            //get all needed handles
            std::array<HalfedgeHandle, 3> boundary_edges;
            HalfedgeHandle heh1 = Base::mesh().next_halfedge_handle(heh);
            HalfedgeHandle heh2 = Base::mesh().prev_halfedge_handle(heh);
            boundary_edges[0] = heh1; boundary_edges[1] = heh2; boundary_edges[2] = heh;

            //calculate E1 and E2
            for (int i = 0; i<3; ++i) {
                heh = boundary_edges[i];        
                VertexHandle vhto = Base::mesh().to_vertex_handle(heh);           
                VertexHandle vhfrom = Base::mesh().from_vertex_handle(heh);
                Vec3d coords0 = vector_cast<Vec3d>(Base::mesh().point(vhto));
                Vec3d coords1 = vector_cast<Vec3d>(Base::mesh().point(vhfrom));
                E1[i] = coords1-coords0;
                E2[i] = coords1.cross(coords0);
                e1 += E1[i];
                e2 += E2[i];
            }
            e3 = e1.cross(e2);

            //equation 7
            constraint = e3*(e1.dot(e1));
            bside = -(e3.dot(e3));
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);
            //equation 8
            constraint = e1.cross(e3);
            bside = 0.0;
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);


            for (int i = 0; i<3; ++i) {
                //create the (e x ) matrices
                e1x[0][1] = -E1[i][2]; e1x[0][2] = E1[i][1]; e1x[1][0] = E1[i][2];
                e1x[1][2] = -E1[i][0]; e1x[2][0] = -E1[i][1]; e1x[2][1] = E1[i][0];
                for (int j = 0; j<3; ++j) {
                    for (int k = 0; k<3; ++k) {
                        Hb[j][k] += e1x[j].dot(e1x[k]);
                    }
                } 
                cb += (E1[i]).cross(E2[i]);
                kb += (E2[i].dot(E2[i]));
            }
            
        }

    //-----------------------------------------------------------------------------------------------------------------------     
    //Volume optimization
        if(Base::mesh().property(LTprops, heh).n != 3) 
            calc_remaining_constraints(heh, Hv, cv);          

    //----------------------------------------------------------------------------------------------------------------------
    //Boundary optimization
        if((Base::mesh().is_boundary(vh0) or Base::mesh().is_boundary(vh1)) and Base::mesh().property(LTprops, heh).n != 3) 
            calc_remaining_constraints(heh, Hb, cb);

    //----------------------------------------------------------------------------------------------------------------------
    //Apply triangle shape opt. if necessary
        if(Base::mesh().property(LTprops, heh).n < 3){
            //insert needed vertices into a set
            for (auto& vertex_handle : {vh0, vh1}) {
                for (typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vertex_handle); vv_it.is_valid(); ++vv_it) {
                    vertex_handles.insert(*vv_it);
                }
            }
            //and erase those, which are not needed
            vertex_handles.erase(vh0);
            vertex_handles.erase(vh1);
            //calculate the Hessian and cs
            for (auto& vertex_handle : vertex_handles){
                tri_shape = vector_cast<Vec3d>(Base::mesh().point(vertex_handle));
                Hs[0][0]++; Hs[1][1]++; Hs[2][2]++;   //add identity matrix
                cs -= tri_shape;
                ks += tri_shape.dot(tri_shape);
            }
            calc_remaining_constraints(heh, Hs, cs);
        }

    //----------------------------------------------------------------------------------------------------------------------
    //Calculate edge collapse error
        if(Base::mesh().property(LTprops, heh).n == 3){
            //get final vertex position
            Matrix3d inv = inverse3x3(Base::mesh().property(LTprops, heh).constraints);
            Vec3d bside = Base::mesh().property(LTprops, heh).b_side;
            //solve system of equations using inverse matrix
            for (int i = 0; i<3; ++i) Base::mesh().property(LTprops, heh).res_vertex_coords[i] = inv[i].dot(bside);
            Vec3d V = Base::mesh().property(LTprops, heh).res_vertex_coords;

            //rescale VertexOptimization variables to match the equation (9)
            for (int i = 0; i<3; ++i) Hv[i] /= 18.0;
            cv /= 18.0;
            kv /= 18.0;
            //rescale BoundaryOptimization variables to match the equation (10)
            for (int i = 0; i<3; ++i) Hb[i] *= 0.5;
            cb *= 0.5;
            kb *= 0.5;
            
            //compute volume and boundary cost
            Vec3d temp1, temp2;
            for (int i = 0; i<3; ++i) {
            temp1[i] = Hv[i].dot(V);
            temp2[i] = Hb[i].dot(V);  
            }
            double fv = 0.5*(V.dot(temp1)) + (cv.dot(V)) + 0.5*kv;  //volume objective function
            double fb = 0.5*(V.dot(temp2)) + (cb.dot(V)) + 0.5*kb;  //area objective function
            //calculate final error
            double length = (vector_cast<Vec3d>(Base::mesh().point(vh1)) - vector_cast<Vec3d>(Base::mesh().point(vh0))).length(); // NAJIT VHODNEJSI FUNKCI EXISTUJE-LI
            double err = lambda*fv +                      //volume opt
                        (1-lambda)*length*length*fb;    //boundary opt


            /*Base::mesh().property(LTprops, rev_heh).res_vertex_coords = Base::mesh().property(LTprops, heh).res_vertex_coords;
            if(A) {
                error.emplace_back(heh);
            }*/
            //Base::mesh().property(LTprops, heh).cost = err; 
            return static_cast<float>(err); 
        } 
    }
    //if 2 or less contraints were found, we won't collapse this edge
    return FLT_MAX;
}

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
preprocess_collapse(const CollapseInfo& _ci)
{
    //move remaining vertex to ideal calculated position
    DefaultTraits::Point ideal_vertex;
    for (int i = 0; i<3; ++i) ideal_vertex[i] = Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords[i];
    Base::mesh().set_point(_ci.v1, ideal_vertex);

    //std::cout<<"Moved vertex to: "<<Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords<<std::endl;
    //std::sort(error.begin(), error.end(), [&](const HalfedgeHandle& a, const HalfedgeHandle& b) { return Base::mesh().property(LTprops, a).e<Base::mesh().property(LTprops, b).e; });
    /*std::ofstream outfile("output.txt", std::ofstream::out | std::ofstream::app);
    outfile << Base::mesh().property(LTprops, _ci.v0v1).cost <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords[0]
    <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords[1]
    <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords[2]
    <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).b_side[0]
    <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).b_side[1]
    <<"\t" <<Base::mesh().property(LTprops, _ci.v0v1).b_side[2];*/
    //"\t"<<Base::mesh().property(LTprops, heh).n<<"\t"<<Base::mesh().property(LTprops, heh).b_side<<"\t"<<Base::mesh().property(LTprops, heh).ncalc<<"\n";
    //outfile<<std::endl;
    //outfile.close();
    //A = false;
}

template<class DecimaterType>
bool
ModLindTurkT<DecimaterType>::
is_alpha_compatible(const HalfedgeHandle& heh, const Vec3d& constraint)
{
    if(Base::mesh().property(LTprops, heh).n == 0){
        return  !(constraint[0] == 0 and constraint[1] == 0 and constraint[2] == 0);
    }
    else if(Base::mesh().property(LTprops, heh).n == 1){
        return (std::pow((Base::mesh().property(LTprops, heh).constraints[0].dot(constraint)), 2)
            < std::pow((Base::mesh().property(LTprops, heh).constraints[0].norm()*constraint.norm()),2)*COSALPHA2);
    }
    else if(Base::mesh().property(LTprops, heh).n == 2){
        Vec3d crossp = Base::mesh().property(LTprops, heh).constraints[0].cross(Base::mesh().property(LTprops, heh).constraints[1]);
        return (std::pow((crossp.dot(constraint)), 2) 
            > std::pow((crossp.norm()*constraint.norm()),2)*SINALPHA2);
    }
    return false;
}

//Adds constraint to the system
template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
add_constraint(const HalfedgeHandle& heh, const Vec3d& constraint, double& right_side)
{
    Base::mesh().property(LTprops, heh).constraints[Base::mesh().property(LTprops, heh).n] = constraint;
    Base::mesh().property(LTprops, heh).b_side[Base::mesh().property(LTprops, heh).n] = right_side;
    Base::mesh().property(LTprops, heh).n++;
}

//Calculates remaining constraints if there are less than 3
template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
calc_remaining_constraints(const HalfedgeHandle& heh, Matrix3d& Hessian, Vec3d& c)
{
    size_t n = Base::mesh().property(LTprops, heh).n;
    size_t N = 3-n;
    //Create identity matrix of a size (3-n, 3)
        std::vector<Vec3d> I(N);
        for (int i = 2-n, j = 2; i>=0; i--, j--) I[i][j] = 1;
    //Create orthogonal matrix Z
        Matrix3d Z;
        Vec3d zeros(0.0);
        Z.fill(zeros);

        for (int i = n-1; i>=0; i--) Z[i] = Base::mesh().property(LTprops, heh).constraints[i];
        if(n == 0) {Z[0][0]++; Z[1][1]++; Z[2][2]++;}                 //If no constraints so far, create a matrix of standard base vectors
        else {
            if (n == 1) {Z[1][0] = Z[0][1]; Z[1][1] = -Z[0][0];}     //Add first orthogonal vector
            Z[2] = Z[0].cross(Z[1]);                                 //Add second orthogonal vector
        }
    //compute remaining constraints and b sides
        std::vector<Vec3d> res(N);
        std::vector<double> resb(N);
        std::vector<Vec3d> temp(N); 
        Matrix3d invZ = inverse3x3(Z);
        for (size_t i = 0; i<3-n; ++i) {
            for (int j = 0; j<3; ++j) {
                temp[i][j] = I[i].dot(invZ[j]);
            }
        }
        for (size_t i = 0; i<3-n; ++i) {
            for (int j = 0; j<3; ++j) {
                res[i][j] = temp[i].dot(Hessian[j]);
            }
        }
        for (size_t i = 0; i<3-n; ++i) resb[i] = -(temp[i].dot(c));
    //add constraints if possible
        for (int i = 2-n; i>=0; --i) if(is_alpha_compatible(heh, res[i])) add_constraint(heh, res[i], resb[i]); 
}

template<class DecimaterType>
Matrix3d
ModLindTurkT<DecimaterType>::inverse3x3(const Matrix3d& m){
    double  a = m[0][0], b = m[1][0], c = m[2][0],
            d = m[0][1], e = m[1][1], f = m[2][1],
            g = m[0][2], h = m[1][2], i = m[2][2];
    double  a00 = e*i-f*h,  a01 = -d*i+f*g, a02 = -g*e+d*h,
            a10 = -b*i+c*h, a11 = a*i-c*g,  a12 = -a*h+b*g,
            a20 = -e*c+f*b, a21 = -a*f+c*d, a22 = a*e-d*b;
    Matrix3d adjugate;
    adjugate[0] = Vec3d(a00, a01, a02);
    adjugate[1] = Vec3d(a10, a11, a12);
    adjugate[2] = Vec3d(a20, a21, a22);
    double det = a*a00+b*a01+c*a02;
    Matrix3d res;
    for (int i = 0; i<3; ++i) res[i] = adjugate[i]*(1.0/det);
    return res;
}


//-----------------------------------------------------------------------------

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
