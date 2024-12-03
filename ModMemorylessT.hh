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


/** \file ModMemorylessT.hh
 */

//=============================================================================
//
//  CLASS ModMemorylessT
//
//=============================================================================

#ifndef OSG_MODMEMORYLESS_HH
#define OSG_MODMEMORYLESS_HH

//== INCLUDES =================================================================

#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <Eigen/Dense>

//== NAMESPACE ================================================================

namespace OpenMesh  {
namespace Decimater {

//== CLASS DEFINITION =========================================================

/** \brief Mesh decimation module computing collapse priority based on 
 *  Memoryless simplification algorithm by Peter Lindstrom and Greg Turk.
 *  (Fast and memory efficient polygonal simplification, 1998) 
 * 
 */
template <class MeshT>
class ModMemorylessT : public ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModMemorylessT, MeshT, Memoryless );

public:

   explicit ModMemorylessT( MeshT &_mesh )
    : Base(_mesh, false)
  {
    // Add needed mesh properties
    Base::mesh().add_property(is_locked);
    Base::mesh().add_property(error_calculated);
    Base::mesh().add_property(n_);
    Base::mesh().add_property(constraints);
    Base::mesh().add_property(rhs);
    Base::mesh().add_property(ideal_vertex_coords);
    Base::mesh().add_property(FProps);
  }
  
  virtual ~ModMemorylessT()
  {
    // Remove properties after simplification is done
    Base::mesh().remove_property(is_locked);
    Base::mesh().remove_property(error_calculated);
    Base::mesh().remove_property(n_);
    Base::mesh().remove_property(constraints);
    Base::mesh().remove_property(rhs);
    Base::mesh().remove_property(ideal_vertex_coords);
    Base::mesh().remove_property(FProps);
  }


public: // inherited

  /** \brief Initalize the module and prepare the mesh for decimation, 
   * possibly lock boundary edges if option is set.
   */
  virtual void initialize(void) override;

  /** \brief Compute ideal vertex position for a halfedge and it's cost
   *  function (error).
   */ 
  virtual float collapse_priority(const CollapseInfo& _ci) override;
	
  /** \brief Move remaining vertex (v1) to the ideal calculated position.
   */
  virtual void preprocess_collapse(const CollapseInfo& _ci) override {
    Base::mesh().set_point(_ci.v1, Base::mesh().property(ideal_vertex_coords, _ci.v0v1));}
  
  /** \brief Recalculate face normals of changed faces (neighbors of v1)
   */
  virtual void postprocess_collapse(const CollapseInfo& _ci) override {
    
    typename Mesh::VertexFaceIter vf_it = Base::mesh().vf_iter(_ci.v1);
    for (; vf_it.is_valid(); ++vf_it) calc_face_normal_and_det(*vf_it);
  }

public: // specific methods

  // Parse simplification parameters string
  void set_opts(std::string opts)
  {
    std::istringstream str(opts);
    std::string cut;

    if (std::getline(str, cut, ',')) {
      // set first parameter (ideal collapse vertex finding mode)
      int cutn = std::stoi(cut);
      if (cutn == 0 or cutn == 1) set_min_mod(cutn);
      else std::cout << "Incorrect value for error minimalization mode (first parameter), going with 0" << "\n"
                     << "Permitted values: 0 and 1" << "\n"
                     << "0 = Original Lindstrom-Turk memoryless simplification,\
                         looking for ideal collapse vertex in 3D space" << "\n"
                     << "1 = Looking for ideal vertex only on line determined by v0 and v1" << std::endl;
    }
    // set second parameter (boundary lock)
    if (std::getline(str, cut, ',')) set_lock(std::stof(cut));
  }

 
private:

  // ------Private functions---------------------------------------------------
  

  /** \brief If the 'lock_boundary_edges' parameter is set to true, the module
   * will preserve the mesh boundary.
   * 
   * \details Check 'initialize' function for details.
   */
  void set_lock(double _lock) { lock_boundary_edges = _lock; }

  /** \brief Checks if constraint is alpha compatible, so we can add it 
   * to the system using function 'add_constraint'.
   * 
   * \param heh The current halfedge
   * \param constr The constraint vector (plane equation vector without
   *  right side)
   * 
   * \return True if the constraint passes all alpha-compatible checks,
   * so it can be added to the system of constraints.
   */
  bool is_alpha_compatible(const HalfedgeHandle& heh, 
                           const Eigen::Vector3d& constr);

  /** \brief Adds a constraint to current halfedge properties and increments
   *  number of constraints for given halfedge.
   */
  void add_constraint(const HalfedgeHandle& heh, 
                      const Eigen::Vector3d& constr, 
                      const double& rhs_);

  /** \brief Calculates remaining constraints if there are less than 3.
   * 
   * \details Using Hessian and c vector which are calculated in volume and
   * boundary preservation sections in the error calculation function, this
   * function tries to find the remaining constraints to fully constraint
   * the ideal collapse vertex to 1 point in the 3D space.
   * 
   * \param heh The current halfedge we are calculating error for
   * \param Hessian The Hessian for given error optimization
   * \param c The vector c for given error optimization
   * 
   * \return Fills the halfedge handle constraints property with found 
   * constraints if they are alpha-compatible.
   */
  void calc_remaining_constraints(const HalfedgeHandle& heh, 
                                  const Eigen::Matrix3d& Hessian, 
                                  const Eigen::Vector3d& c);

  /** \brief Calculates face normal and determinant of a matrix
   *  determined by 3 face vertices (columns).
   */
  void calc_face_normal_and_det(const FaceHandle& fh);


  /** \brief Function for casting vertex coordinates to Eigen vector
   */
  inline Eigen::Vector3d eigenvec_cast(const VertexHandle& v) {
    DefaultTraits::Point p = Base::mesh().point(v);
    return Eigen::Vector3d(p[0], p[1], p[2]);
  }



  // ------Parameters of the decimating module---------------------------------


  /** Parameter for locking all boundary and "semi-boundary" edge (edges with 
   * only one boundary vertex) to preserve the mesh boundary. If we don't
   * collapse these edges, the boundary will stay the same.
   */
  bool lock_boundary_edges = false;  

  // constraint compatibility parameters
  double  SINALPHA = std::sin(0.01745329251),
          COSALPHA = std::cos(0.01745329251);

  // Defines, how will the module find ideal collapse vertex
  // SPACE = original LT, searches the whole 3D space
  // LINE = search restricted to only one line (edge v0v1) (= 2 constraints already given)
  enum decimation_mode { SPACE = 0, LINE = 1 };
  decimation_mode mod = SPACE;

  /** \brief Sets the mode for finding ideal collapse vertex
   */
  void set_min_mod(const int& mod_) {
    switch(mod_){
      case 0: mod = SPACE;break;
      case 1: mod = LINE; break;
      default:mod = SPACE;break;
    }
  }



  // ------Properties of each halfedge-----------------------------------------


  /** If the halfedge is locked, the error will be automatically FLT_MAX.
   *  - This property is activated with the parameter "lock_boundary_edges"
   * -> if set to true, it will set this property to true to every boundary
   * and "semi-boudary" edges (more concretely their respective halfedges)
   */
  HPropHandleT<bool> is_locked;

  /** If the error is already calculated for the opposite halfedge,
   * no need to calculate it, so let's set it to FLT_MAX
   * (the error calculating function checks if the opposite halfedge has
   * it's error calculated: if yes, it sets the error for the current halfedge
   * to FLT_MAX)
   */
  HPropHandleT<bool> error_calculated;

  /** Number of valid constraints acquired for given halfedge during the error 
   * calculation.
   * (number of rows in matrix A)
   */
  HPropHandleT<size_t> n_;

  /** Storage for all valid constraints acquired for given halfedge during the 
   * error calculation. (= matrix A)
   */
  HPropHandleT<Eigen::Matrix3d> constraints;

  // The right side of equation Av=b
  HPropHandleT<Eigen::Vector3d> rhs;

  /** Ideal vertex position after collapse
   * Solution to equation Av=b solved using: v = A^(-1) * b
   */
  HPropHandleT<DefaultTraits::Point> ideal_vertex_coords;

  // Storage for basic face properties to save repeating calculations
  struct FaceProps {

    double det;                         // determinant of vertices of the face
    double det2;                        // determinant^2
    Eigen::Vector3d face_normal;        // classic face normal
    Eigen::Vector3d det_dot_normal;     // determinant*normal
    Eigen::Matrix3d face_normal_matrix; // normal*normal^T
  };
  // Declaring the struct as a face property
  FPropHandleT<FaceProps> FProps;
};

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_DECIMATER_MODMEMORYLESS_CC)
#define OSG_MODMEMORYLESS_TEMPLATES
#include "ModMemorylessT_impl.hh"
#endif
//=============================================================================
#endif // OSG_MODMEMORYLESS_HH defined
//=============================================================================
