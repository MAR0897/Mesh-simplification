/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2025, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */


//=============================================================================
//
//  CLASS ModLindTurk - IMPLEMENTATION
//
//=============================================================================

#define OPENMESH_DECIMATER_MODLINDTURK_CC

//== INCLUDES =================================================================

#include <OpenMesh/Tools/Decimater/ModLindTurkT.hh>

//== NAMESPACE ================================================================

namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Decimater { // BEGIN_NS_DECIMATER

//== IMPLEMENTATION ===========================================================

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
initialize()
{
    if (!is_locked.is_valid()) Base::mesh().add_property(is_locked);
    if (!error_calculated.is_valid()) Base::mesh().add_property(error_calculated);
    if (!n_.is_valid()) Base::mesh().add_property(n_);
    if (!constraints.is_valid()) Base::mesh().add_property(constraints);
    if (!rhs.is_valid()) Base::mesh().add_property(rhs);
    if (!ideal_vertex_coords.is_valid()) Base::mesh().add_property(ideal_vertex_coords);
    if (!FProps.is_valid()) Base::mesh().add_property(FProps);
    

    typename Mesh::HalfedgeIter he_it = Base::mesh().halfedges_begin(),
                                he_end = Base::mesh().halfedges_end();

    // initialize bool halfedge properties 
    for (; he_it != he_end; ++he_it) {
        Base::mesh().property(error_calculated, *he_it) = false;
        Base::mesh().property(is_locked, *he_it) = false;
    }

    //Lock all boundary edges if option is set, lock all boundary
    // and "semi-boundary" edges, so that the mesh boundary stays the same
    he_it = Base::mesh().halfedges_begin();
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
                    Base::mesh().property(is_locked, *voh_it1) = true;
                for (; voh_it2.is_valid(); ++voh_it2) 
                    Base::mesh().property(is_locked, *voh_it2) = true;
                for (; vih_it1.is_valid(); ++vih_it1) 
                    Base::mesh().property(is_locked, *vih_it1) = true;
                for (; vih_it2.is_valid(); ++vih_it2) 
                    Base::mesh().property(is_locked, *vih_it2) = true;
            }
        }
    }

    // calculate face normals and store them
    typename Mesh::FaceIter f_it = Base::mesh().faces_begin(),
                            f_end = Base::mesh().faces_end();

    for (; f_it != f_end; ++f_it) calc_face_normal_and_det(*f_it);
}

//=============================================================================

template<class DecimaterType>
float
ModLindTurkT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci) 
{
    // only take not locked halfedge, and halfedges, whose opposite halfedge has not calculated error yet
    if(!Base::mesh().property(is_locked, _ci.v0v1) and !Base::mesh().property(error_calculated, _ci.v1v0)){

    //-------------------------------------------------------------------------
    // 0. Initializing variables
    //-------------------------------------------------------------------------
        HalfedgeHandle heh = _ci.v0v1;    //halfedge we are currently calculating error for

        bool v0_is_boundary = Base::mesh().is_boundary(_ci.v0);
        bool v1_is_boundary = Base::mesh().is_boundary(_ci.v1);
        bool is_boundary = (v0_is_boundary or v1_is_boundary); // true for boundary and semi-boundary edges

        //set variables to zero
        Base::mesh().property(n_, heh) = 0;                      
        Base::mesh().property(constraints, heh).setZero();
        Base::mesh().property(rhs, heh).setZero();           

        Eigen::Vector3d constraint; constraint.setZero();   //storage for new constraint
        double bside = 0.0;                                 //storage for new right side number

        std::set<FaceHandle> face_handles;
           
        // Hessian, c vector and k for final error calculation
        // (and possibly for calculating remaining constraints)
        Eigen::Matrix3d Hv, Hb; Hv.setZero();
        Eigen::Vector3d cv, cb; cv.setZero();
        double kv = 0.0, kb = 0.0;

    //-------------------------------------------------------------------------
    // 0,5. Adding first two constraints if decimation mode is LINE
    //-------------------------------------------------------------------------
        if (mod == LINE) {

            Eigen::Vector3d a1, a2; a1.setZero(); a2.setZero();
            double b1 = 0.0, b2 = 0.0;

            double dx = _ci.p1[0] - _ci.p0[0];
            double dy = _ci.p1[1] - _ci.p0[1];
            double dz = _ci.p1[2] - _ci.p0[2];
            Eigen::Vector3d A = eigenvec_cast(_ci.v0);


            // create perpendicular vectors to halfedge vector
            if (dx == 0.0) {
                a1 = {1.0,0.0,0.0};
                if (dy == 0.0) a2 = {0.0, 1.0, 0.0};
                else if (dz == 0.0) a2 = {0.0, 0.0, 1.0};
                else a2 = {0.0, -dz, dy};
            }
            else if (dy == 0.0){
                a1 = {0.0,1.0,0.0};
                if (dz == 0.0) a2 = {0.0, 0.0, 1.0};
                else a2 = {-dz, 0.0, dx};
            }
            else if (dz == 0.0) { 
                a1 = {0.0,0.0,1.0}; 
                a2 = {-dy, dx, 0.0}; 
            } 
            else { 
                a1 = {-dz, 0.0, dx}; 
                a2 = {-dy, dx, 0.0};
            }

            b1 = a1.dot(A);
            b2 = a2.dot(A);

            if(is_alpha_compatible(heh, a1)) add_constraint(heh, a1, b1);
            if(is_alpha_compatible(heh, a2)) add_constraint(heh, a2, b2);
            else {  // if dx is too small, the two constraints are not alpha compatible, 
                    // so we need to create one without dx coordinate
                a2 = Eigen::Vector3d{0.0, -dz, dy};
                b2 = a2.dot(A);
                if(is_alpha_compatible(heh, a2)) add_constraint(heh, a2, b2);
            }
        }

    //-------------------------------------------------------------------------
    // 1. Volume preservation (+volume optimization)
    //-------------------------------------------------------------------------
        // plug the faces we need into a set
        for (auto& vh : {_ci.v0, _ci.v1}) {
            typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(vh);
            for (; vf_it.is_valid(); ++vf_it) face_handles.insert(*vf_it);
        }

        // calc constraint and Vertex Optimization variables
        for (const auto& ff : face_handles) {
            
            constraint += Base::mesh().property(FProps, ff).face_normal;
            bside += Base::mesh().property(FProps, ff).det;

            Hv += Base::mesh().property(FProps, ff).face_normal_matrix;
            cv -= Base::mesh().property(FProps, ff).det_dot_normal;
            kv += Base::mesh().property(FProps, ff).det2;
        }
        //rescale Vertex Optimization variables to match the equation (9)
        Hv /= 18.0; cv /= 18.0; kv /= 36.0;  
       
        if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);
    //-------------------------------------------------------------------------
    // 2. Boundary preservation (+boundary optimization)
    //-------------------------------------------------------------------------
        if(is_boundary and (Base::mesh().property(n_, heh) < 3)){

            // if edge is boundary, there will be 3 edges needed for constraints
            // calculation (Figure 3), if not, there will be only 2
            size_t N;
            if (Base::mesh().is_boundary(heh)) N = 3; else N = 2;
            std::vector<HalfedgeHandle> boundary_edges;

            Hb.setZero(); cb.setZero(); kb = 0.0;

            if (N == 3) {
                HalfedgeHandle heh1 = Base::mesh().next_halfedge_handle(heh);
                HalfedgeHandle heh2 = Base::mesh().prev_halfedge_handle(heh);
                for (const auto& hh : {heh1, heh2, heh}) boundary_edges.emplace_back(hh); 
            }
            else {
                //getting the two boundary edges
                typename Mesh::VertexOHalfedgeIter voh_it;
                typename Mesh::VertexIHalfedgeIter vih_it;
                //if v0 is boundary, check its outgoing halfedges, if not, check v1 outgoing halfedges
                if (v0_is_boundary) {
                    voh_it = Base::mesh().voh_iter(_ci.v0); 
                    vih_it = Base::mesh().vih_iter(_ci.v0);
                } 
                else {
                    voh_it = Base::mesh().voh_iter(_ci.v1); 
                    vih_it = Base::mesh().vih_iter(_ci.v1);
                }
                for (; voh_it.is_valid(); ++voh_it) 
                    if (Base::mesh().is_boundary(*voh_it)) 
                        boundary_edges.emplace_back(*voh_it);
                for (; vih_it.is_valid(); ++vih_it) 
                    if (Base::mesh().is_boundary(*vih_it)) 
                        boundary_edges.emplace_back(*vih_it);
            }

            Eigen::Matrix3d E1, E2, e1x;  E1.setZero(); E2.setZero(); e1x.setZero();               
            Eigen::Vector3d e1, e2, e3; e1.setZero(); e2.setZero(); e3.setZero();                     

            //calculate E1 and E2 for every edge (that is, for 2 or 3 edges)
            for (size_t i = 0; i<N; ++i) {
                HalfedgeHandle hh = boundary_edges[i];        
                Eigen::Vector3d coords0 = eigenvec_cast(Base::mesh().from_vertex_handle(hh));
                Eigen::Vector3d coords1 = eigenvec_cast(Base::mesh().to_vertex_handle(hh));
                E1.row(i) = coords1-coords0;
                E2.row(i) = coords1.cross(coords0);
                e1 += E1.row(i);
                e2 += E2.row(i);
            }
            e3 = e1.cross(e2);
            
            // equation 7
            constraint = e3*(e1.transpose()*e1);   bside = -(e3.transpose()*e3).value();
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);
            // equation 8
            constraint = e1.cross(e3);             bside = 0.0;
            if(is_alpha_compatible(heh, constraint)) add_constraint(heh, constraint, bside);    

            //calculate hessian, c and k from values computed at the boundary preservation section
            for (size_t i = 0; i<N; ++i) {
                //create the (e x ) matrices
                e1x(0,1) = -E1(i,2); e1x(0,2) = E1(i,1); e1x(1,0) = E1(i,2);
                e1x(1,2) = -E1(i,0); e1x(2,0) = -E1(i,1); e1x(2,1) = E1(i,0);
                Hb += e1x*e1x.transpose();
                cb += (E1.row(i)).cross(E2.row(i));
                kb += (E2.row(i)).dot(E2.row(i));
            }

            //rescale Boundary Optimization variables to match the equation (10)
            Hb *= 0.5; cb *= 0.5; kb *= 0.25;     
        }
    //-------------------------------------------------------------------------
    // 3. Volume optimization
    //-------------------------------------------------------------------------
        if(Base::mesh().property(n_, heh) < 3) 
            calc_remaining_constraints(heh, Hv, cv);
    //-------------------------------------------------------------------------
    // 4. Boundary optimization
    //-------------------------------------------------------------------------
        if(is_boundary and Base::mesh().property(n_, heh) < 3) 
            calc_remaining_constraints(heh, Hb, cb);
    //-------------------------------------------------------------------------
    // 5. Triangle shape optimization
    //-------------------------------------------------------------------------
        if(Base::mesh().property(n_, heh) < 3){
            //std::cout<<III++<<"\t";
            std::set<VertexHandle> vertex_handles;
            Eigen::Matrix3d Hs; Hs.setZero();
            Eigen::Vector3d cs; cs.setZero();
            double ks = 0.0;

            //insert needed vertices into a set
            for (auto& vh : {_ci.v0, _ci.v1}) {
                typename Mesh::VertexVertexIter vv_it = Base::mesh().vv_iter(vh);
                for (; vv_it.is_valid(); ++vv_it) vertex_handles.insert(*vv_it);
            }
            //and erase those, which are not needed
            vertex_handles.erase(_ci.v0);
            vertex_handles.erase(_ci.v1);

            //calculate the Hessian and cs
            for (auto& vv : vertex_handles){
                Eigen::Vector3d tri_shape = eigenvec_cast(vv);
                Hs(0,0)+=2; Hs(1,1)+=2; Hs(2,2)+=2;   //add identity matrix
                cs -= 2*tri_shape;
                ks += 2*(tri_shape.dot(tri_shape));
            }
            calc_remaining_constraints(heh, Hs, cs);
        }
    //-------------------------------------------------------------------------
    // 6. Calculate edge collapse error
    //-------------------------------------------------------------------------
        if(Base::mesh().property(n_, heh) == 3){

            // get final vertex position (solve system of equations using inverse matrix) Av = b
            Eigen::Vector3d V = Base::mesh().property(constraints, heh).inverse()*
                                Base::mesh().property(rhs, heh);

            // store the vertex position as Point
            for (int i = 0; i<3; ++i) Base::mesh().property(ideal_vertex_coords, heh)[i] = V[i];
        
            // compute volume and boundary cost
            double fv, fb = 0.0;   fv = (0.5*(V.transpose()*(Hv*V)) + (cv.transpose()*V)).value() + kv;   //volume objective function
            if (is_boundary)       fb = (0.5*(V.transpose()*(Hb*V)) + (cb.transpose()*V)).value() + kb;   //area objective function
 
            // final error E = lambda*fv + (1-lambda)*L^2*fb
            double err = 0.0;
            if (fb == 0.0) err = 0.5*fv;
            else {
                double length = (eigenvec_cast(_ci.v1)-eigenvec_cast(_ci.v0)).norm();
                err = 0.5*(fv + length*length*fb);
            }
        
            Base::mesh().property(error_calculated, _ci.v0v1) = true;

            return static_cast<float>(err); 
        } 

    }

    // if 2 or less contraints were found, we won't collapse this edge
    // (or if the halfedge was locked, or its opposite has already calculated error)
    return FLT_MAX;
}

//=============================================================================

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
calc_face_normal_and_det(const FaceHandle& fh) 
{
    Eigen::Matrix3d fv_coords;
    typename Mesh::FaceVertexIter fv_it = Base::mesh().fv_iter(fh);
    for (size_t i = 0; fv_it.is_valid(); ++fv_it, ++i)
        fv_coords.col(i) = eigenvec_cast(*fv_it);

    Eigen::Vector3d AB = fv_coords.col(1)-fv_coords.col(0);
    Eigen::Vector3d AC = fv_coords.col(2)-fv_coords.col(0);
    Eigen::Vector3d normal = AB.cross(AC);
    double determinant = fv_coords.col(0).dot(normal);
    
    Base::mesh().property(FProps, fh).face_normal = normal;
    Base::mesh().property(FProps, fh).det = determinant;
    Base::mesh().property(FProps, fh).det2 = determinant*determinant;
    Base::mesh().property(FProps, fh).det_dot_normal = determinant*normal;
    Base::mesh().property(FProps, fh).face_normal_matrix = normal*normal.transpose();
}

//=============================================================================

template<class DecimaterType>
bool
ModLindTurkT<DecimaterType>::
is_alpha_compatible(const HalfedgeHandle& heh, const Eigen::Vector3d& constr)
{   
    // a1 != null vector
    if(Base::mesh().property(n_, heh) == 0)
        return  !(constr(0) == 0.0 and constr(1) == 0.0 and constr(2) == 0.0);

    // ((a1^T)*a2)^2 < (||a1||*||a2||*cos(alpha))^2
    else if(Base::mesh().property(n_, heh) == 1){
        Eigen::Vector3d a1 = Base::mesh().property(constraints, heh).row(0);
        return (std::abs(a1.dot(constr)) < std::abs((a1.norm()*constr.norm())*COSALPHA));
    }
    // ((a1 x a2)^T * a3)^2 > (||a1 x a2||*||a3||*sin(alpha))^2
    else if(Base::mesh().property(n_, heh) == 2){
        Eigen::Vector3d crossp =  Base::mesh().property(constraints, heh).row(0).cross(Base::mesh().property(constraints, heh).row(1));
        return (std::abs((crossp.dot(constr))) > std::abs((crossp.norm()*constr.norm())*SINALPHA));
    }
    return false;
}

//=============================================================================

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
add_constraint( const HalfedgeHandle& heh, 
                const Eigen::Vector3d& constr, 
                const double& rhs_)
{
    Base::mesh().property(constraints, heh).row(Base::mesh().property(n_, heh)) = constr;
    Base::mesh().property(rhs, heh)[Base::mesh().property(n_, heh)] = rhs_;
    Base::mesh().property(n_, heh)++;
}

//=============================================================================

template<class DecimaterType>
void
ModLindTurkT<DecimaterType>::
calc_remaining_constraints( const HalfedgeHandle& heh,
                            const Eigen::Matrix3d& Hessian, 
                            const Eigen::Vector3d& c)
{
    size_t n = Base::mesh().property(n_, heh);
 
    //Create orthogonal matrix Z
    Eigen::Matrix3d Z = Base::mesh().property(constraints, heh).transpose();
    if(n == 0) Z = Eigen::MatrixXd::Identity(3,3);        //If no constraints so far, create a matrix of standard base vectors
    else {
        if (n == 1) {Z(0,1) = Z(1,0); Z(1,1) = -Z(0,0);}  //Add first orthogonal vector
        Z.col(2) = Z.col(0).cross(Z.col(1));              //Add second orthogonal vector
    }
  
    //compute remaining constraints and b sides        
    auto temp = (Z.inverse()).block(n, 0, 3-n, 3);
    auto constraints = temp*Hessian;
    auto bsides = -temp*c;
    
    //add constraints if possible
    for (size_t i = 0; i<=2-n; ++i)
        if(is_alpha_compatible(heh, constraints.row(i)))       
            add_constraint(heh, constraints.row(i), bsides(i));
}

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
