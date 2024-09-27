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
//  CLASS ModQuadricT
//
//=============================================================================

#ifndef OSG_MODQUADRIC_HH
#define OSG_MODQUADRIC_HH


//== INCLUDES =================================================================

#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Geometry/QuadricT.hh>
#include <OpenMesh/Eigen/Dense>


//== NAMESPACE ================================================================

namespace OpenMesh  {
namespace Decimater {


//== CLASS DEFINITION =========================================================


/** \brief Mesh decimation module computing collapse priority based on error quadrics.
 *
 *  This module can be used as a binary and non-binary module.
 */
template <class MeshT>
class ModQuadricT : public ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModQuadricT, MeshT, Quadric );

public:

  /** Constructor
   *  \internal
   */
  explicit ModQuadricT( MeshT &_mesh )
    : Base(_mesh, false)
  {
    unset_max_err();
    Base::mesh().add_property( quadrics_ );
    Base::mesh().add_property(res_vertex_coords);
  }


  /// Destructor
  virtual ~ModQuadricT()
  {
    Base::mesh().remove_property(quadrics_);
    Base::mesh().remove_property(res_vertex_coords);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation.
  virtual void initialize(void) override;

  /** Compute collapse priority based on error quadrics.
   *
   *  \see ModBaseT::collapse_priority() for return values
   *  \see set_max_err()
   */
  virtual float collapse_priority(const CollapseInfo& _ci) override
  {
    using namespace OpenMesh;

    typedef Geometry::QuadricT<double> Q;

    Q q = Base::mesh().property(quadrics_, _ci.v0);
    q += Base::mesh().property(quadrics_, _ci.v1);


    //choose a perfect spot for vertex v (calculate ideal coords)
    Eigen::Matrix4d V;
    Eigen::Vector4d v_coords;
    V.setZero();
    v_coords.setZero();
    V(0,0) = q.a();
    V(0,1) = q.b();
    V(0,2) = q.c();
    V(0,3) = q.d();
    V(1,1) = q.e();
    V(1,2) = q.f();
    V(1,3) = q.g();
    V(2,2) = q.h();
    V(2,3) = q.i();
    V(3,3) = q.j();
    V(3,0) = 0; V(3,1) = 0; V(3,2) = 0; V(3,3) = 1;
    v_coords = V.inverse()*(Eigen::Vector4d{0.0, 0.0, 0.0, 1.0});
    for (size_t j = 0; j<3; ++j) Base::mesh().property(res_vertex_coords, _ci.v0v1)[j] = v_coords[j];



    double err = v_coords.transpose()*(V*v_coords); //q(_ci.p1);

    //min_ = std::min(err, min_);
    //max_ = std::max(err, max_);

    //double err = q( p );

    return float( (err < max_err_) ? err : float( Base::ILLEGAL_COLLAPSE ) );
  }

  
  virtual void preprocess_collapse(const CollapseInfo& _ci)
  {
      //move remaining vertex to ideal calculated position
      DefaultTraits::Point ideal_vertex;
      for (int i = 0; i<3; ++i) ideal_vertex[i] = Base::mesh().property(res_vertex_coords, _ci.v0v1)[i];
      Base::mesh().set_point(_ci.v1, ideal_vertex);
  }
  

  /// Post-process halfedge collapse (accumulate quadrics)
  virtual void postprocess_collapse(const CollapseInfo& _ci) override
  {
    Base::mesh().property(quadrics_, _ci.v1) +=
      Base::mesh().property(quadrics_, _ci.v0);
  }

  /// set the percentage of maximum quadric error
  void set_error_tolerance_factor(double _factor) override;



public: // specific methods

  /** Set maximum quadric error constraint and enable binary mode.
   *  \param _err    Maximum error allowed
   *  \param _binary Let the module work in non-binary mode in spite of the
   *                 enabled constraint.
   *  \see unset_max_err()
   */
  void set_max_err(double _err, bool _binary=true)
  {
    max_err_ = _err;
    Base::set_binary(_binary);
  }

  /// Unset maximum quadric error constraint and restore non-binary mode.
  /// \see set_max_err()
  void unset_max_err(void)
  {
    max_err_ = DBL_MAX;
    Base::set_binary(false);
  }

  /// Return value of max. allowed error.
  double max_err() const { return max_err_; }


private:

  // maximum quadric error
  double max_err_;

  // this vertex property stores a quadric for each vertex
  VPropHandleT< Geometry::QuadricT<double> >  quadrics_;

  HPropHandleT<Vec3d> res_vertex_coords;
};

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_DECIMATER_MODQUADRIC_CC)
#define OSG_MODQUADRIC_TEMPLATES
#include "ModQuadricT_impl.hh"
#endif
//=============================================================================
#endif // OSG_MODQUADRIC_HH defined
//=============================================================================

