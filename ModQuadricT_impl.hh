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



/** \file ModQuadricT_impl.hh
    Bodies of template member function.
 */

//=============================================================================
//
//  CLASS ModQuadric - IMPLEMENTATION
//
//=============================================================================

#define OPENMESH_DECIMATER_MODQUADRIC_CC

//== INCLUDES =================================================================

#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>

//== NAMESPACE ===============================================================

namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Decimater { // BEGIN_NS_DECIMATER

//== IMPLEMENTATION ==========================================================


template<class DecimaterType>
void
ModQuadricT<DecimaterType>::
initialize()
{

  // print simplification parameters
  std::cout<< "Garland-Heckbert quadric simplification options"<< "\n"
           << " - " << "vertex selection mode = "<< ((mod == SPACE) ? "SPACE" :
                                                (mod == LINE) ? "LINE" :
                                                (mod == POINTS) ? "POINTS" :
                                                "DEFAULT_OpenMesh") << "\n"
           << " - " << "boundary lock = "<< ((lock_boundary_edges) ? 
                                           "True (mesh boundary wont move)" : "False") <<"\n"
           << " - " << "max error = "<< max_err_ << "\n"
           << std::endl;
  

  using Geometry::Quadricd;
  // alloc quadrics
  if (!quadrics_.is_valid()){
    Base::mesh().add_property( quadrics_ );
    Base::mesh().add_property( ideal_vertex_coords );
    Base::mesh().add_property( is_locked );
    Base::mesh().add_property( error_calculated );}


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

  // clear quadrics
  typename Mesh::VertexIter  v_it  = Base::mesh().vertices_begin(),
                             v_end = Base::mesh().vertices_end();

  for (; v_it != v_end; ++v_it)
    Base::mesh().property(quadrics_, *v_it).clear();

  // calc (normal weighted) quadric
  typename Mesh::FaceIter          f_it  = Base::mesh().faces_begin(),
                                   f_end = Base::mesh().faces_end();

  typename Mesh::FaceVertexIter    fv_it;
  typename Mesh::VertexHandle      vh0, vh1, vh2;
  typedef Vec3d                    Vec3;

  for (; f_it != f_end; ++f_it)
  {
    fv_it = Base::mesh().fv_iter(*f_it);
    vh0 = *fv_it;  ++fv_it;
    vh1 = *fv_it;  ++fv_it;
    vh2 = *fv_it;

    Vec3 v0, v1, v2;
    {
      using namespace OpenMesh;

      v0 = vector_cast<Vec3>(Base::mesh().point(vh0));
      v1 = vector_cast<Vec3>(Base::mesh().point(vh1));
      v2 = vector_cast<Vec3>(Base::mesh().point(vh2));
    }

    Vec3 n = (v1-v0) % (v2-v0);
    double area = n.norm();
    if (area > FLT_MIN)
    {
      n /= area;
      area *= 0.5;
    }

    const double a = n[0];
    const double b = n[1];
    const double c = n[2];
    const double d = -(vector_cast<Vec3>(Base::mesh().point(vh0))|n);

    Quadricd q(a, b, c, d);
    q *= area;

    Base::mesh().property(quadrics_, vh0) += q;
    Base::mesh().property(quadrics_, vh1) += q;
    Base::mesh().property(quadrics_, vh2) += q;
  }
}


//-----------------------------------------------------------------------------

template<class DecimaterType>
float
ModQuadricT<DecimaterType>::
collapse_priority(const CollapseInfo& _ci)
  {
    using namespace OpenMesh;

    //process only not-locked edges and those edges, which dont have computed error yet
    if (!Base::mesh().property(is_locked, _ci.v0v1) and !Base::mesh().property(error_calculated, _ci.v1v0)) {

      typedef Geometry::QuadricT<double> Q;

      // Q = Q1 + Q2
      Q q = Base::mesh().property(quadrics_, _ci.v0);
      q += Base::mesh().property(quadrics_, _ci.v1);
       // final edge error
      double err;
      
      //=======================================================================

     
      if (mod == SPACE) { //choose a ideal collapse vertex coords in 3D space
      
        Eigen::Matrix4d V;        V.setZero();
        Eigen::Vector4d v_coords; v_coords.setZero();
        
        // define matrix V
        V(0,0) = q.a(); V(0,1) = q.b(); V(0,2) = q.c(); V(0,3) = q.d();
        V(1,0) = q.b(); V(1,1) = q.e(); V(1,2) = q.f(); V(1,3) = q.g();
        V(2,0) = q.c(); V(2,1) = q.f(); V(2,2) = q.h(); V(2,3) = q.i();
                                                        V(3,3) = 1.0;
        // calculate ideal coords
        v_coords = V.inverse()*(Eigen::Vector4d{0.0, 0.0, 0.0, 1.0}); 
        // assign ideal coords to a mesh property
        for (size_t j = 0; j<3; ++j) 
          Base::mesh().property(ideal_vertex_coords, _ci.v0v1)[j] = v_coords[j];
        // compute collapse cost (error)
        err = compute_error(q, std::move(v_coords[0]), std::move(v_coords[1]), std::move(v_coords[2]));
      }

      //-----------------------------------------------------------------------

      else if (mod == LINE) { // compute the current edge parametric function and minimize for t
               
        // differences between vertices
        double dx = _ci.p1[0] - _ci.p0[0];
        double dy = _ci.p1[1] - _ci.p0[1];
        double dz = _ci.p1[2] - _ci.p0[2];
        // v0 coords
        double vx = _ci.p0[0];
        double vy = _ci.p0[1];
        double vz = _ci.p0[2];

        // find the minimum on line defined by the current edge
        double t = - (dx*vx*q.a() + dy*vy*q.e() + dz*vz*q.h() + 
                  dx*q.d() + dy*q.g() + dz*q.i() + 
                  q.b()*(dx*vy+dy*vx) + q.c()*(dx*vz+dz*vx) + q.f()*(dz*vy+dy*vz))
                  /
                  (dx*dx*q.a() + dy*dy*q.e() + dz*dz*q.h() +
                  2*q.b()*dx*dy + 2*q.c()*dx*dz + 2*q.f()*dz*dy );

        Eigen::Vector4d v_coords; v_coords.setZero();
        // calculate ideal collapse vertex coords using parametric equation of the line
        v_coords[0] = dx*t + vx;
        v_coords[1] = dy*t + vy;
        v_coords[2] = dz*t + vz;
        // assign ideal coords to a mesh property
        for (size_t j = 0; j<3; ++j) 
          Base::mesh().property(ideal_vertex_coords, _ci.v0v1)[j] = v_coords[j];
        // compute collapse cost (error)
        err = compute_error(q, std::move(v_coords[0]), std::move(v_coords[1]), std::move(v_coords[2]));        
      }

      //-----------------------------------------------------------------------

      else if (mod == POINTS) { // compute errors only in v0, v1 and midpoint of the current edge

        // find the midpoint of currect halfedge
        double midp[3] = {(_ci.p0[0]+_ci.p1[0])/2.0, 
                          (_ci.p0[1]+_ci.p1[1])/2.0, 
                          (_ci.p0[2]+_ci.p1[2])/2.0};
        double midpp[3];
        for (size_t j = 0; j<3; ++j) midpp[j] = midp[j];

        // compute error in v0, v1 and the midpoint
        double err1 = compute_error(q, _ci.p0[0], _ci.p0[1], _ci.p0[2]);
        double err2 = compute_error(q, _ci.p1[0], _ci.p1[1], _ci.p1[2]);
        double err3 = compute_error(q, std::move(midpp[0]), std::move(midpp[1]), std::move(midpp[2]));

        // compare the errors
        if (err1 < err2 and err1 < err3) {
          err = err1;
          Base::mesh().property(ideal_vertex_coords, _ci.v0v1) = _ci.p0;
        }
        else if (err2 < err3) {
          err = err2;
          Base::mesh().property(ideal_vertex_coords, _ci.v0v1) = _ci.p1;
        }
        else {
          err = err3;
          for (size_t j = 0; j<3; ++j) 
            Base::mesh().property(ideal_vertex_coords, _ci.v0v1)[j] = midp[j];
        } 
      }

      //-----------------------------------------------------------------------

      else { // use original OpenMesh implementation, which calculates error only for v0
        
        err = q(_ci.p1);
        //min_ = std::min(err, min_);
        //max_ = std::max(err, max_);
        //double err = q( p );
      }

      //=======================================================================

      // error is calculated, so update property and return error
      if (mod != 0) Base::mesh().property(error_calculated, _ci.v0v1) = true;
      return float( (err < max_err_) ? err : float( Base::ILLEGAL_COLLAPSE ) );
    }

    return FLT_MAX;
  }


//-----------------------------------------------------------------------------

template<class DecimaterType>
double
ModQuadricT<DecimaterType>::
compute_error(Geometry::QuadricT<double>& q, double&& x, double&& y, double&& z) {

  return  x * (x*q.a() + 2*(y*q.b()+q.d())) + 
          y * (y*q.e() + 2*(z*q.f()+q.g())) + 
          z * (z*q.h() + 2*(x*q.c()+q.i())) + q.j();   
}


//-----------------------------------------------------------------------------

template<class MeshT>
void ModQuadricT<MeshT>::set_error_tolerance_factor(double _factor) {
  if (this->is_binary()) {
    if (_factor >= 0.0 && _factor <= 1.0) {
      // the smaller the factor, the smaller max_err_ gets
      // thus creating a stricter constraint
      // division by error_tolerance_factor_ is for normalization
      double new_max_err = max_err_ * _factor / this->error_tolerance_factor_;
      set_max_err(new_max_err);
      this->error_tolerance_factor_ = _factor;

      initialize();
    }
  }
}

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
