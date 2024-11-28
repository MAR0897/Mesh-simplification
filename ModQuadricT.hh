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
#include <Eigen/Dense>


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
    Base::mesh().add_property(ideal_vertex_coords);
    Base::mesh().add_property(is_locked);
    Base::mesh().add_property(error_calculated);
  }

  /// Destructor
  virtual ~ModQuadricT()
  {
    Base::mesh().remove_property(quadrics_);
    Base::mesh().remove_property(ideal_vertex_coords);
    Base::mesh().remove_property(is_locked);
    Base::mesh().remove_property(error_calculated);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation.
  virtual void initialize(void) override;

  /** Compute collapse priority based on error quadrics.
   *
   *  \see ModBaseT::collapse_priority() for return values
   *  \see set_max_err()
   */
  virtual float collapse_priority(const CollapseInfo& _ci) override;

  double compute_error(Geometry::QuadricT<double>& q, double&& x, double&& y, double&& z);
  
  //move remaining vertex to ideal calculated position
  virtual void preprocess_collapse(const CollapseInfo& _ci)
  {
      Base::mesh().set_point(_ci.v1, Base::mesh().property(ideal_vertex_coords, _ci.v0v1));
  }
  
  /// Post-process halfedge collapse (accumulate quadrics)
  virtual void postprocess_collapse(const CollapseInfo& _ci) override
  {
      Base::mesh().property(quadrics_, _ci.v1) +=
      Base::mesh().property(quadrics_, _ci.v0);

      //restart the error calculation property
      /*if (min_mod != 0) {
        typename Mesh::HalfedgeIter he_it = Base::mesh().halfedges_begin(),
                                    he_end = Base::mesh().halfedges_end();
        for (; he_it != he_end; ++he_it) Base::mesh().property(error_calculated, *he_it) = false;
      }*/
  }

  /// set the percentage of maximum quadric error
  void set_error_tolerance_factor(double _factor) override;

private:

  // Defines how is the ideal collapse vertex being calculated
  enum decimation_mode {
    SPACE = 3,      // 3 = original GH (choosing the ideal vertex anywhere from the 3D space)
    LINE = 2,       // 2 = choosing ideal vertex on line determined by v0 and v1
    POINTS = 1,     // 1 = choosing ideal vertex at start/end/mid points only;
    DEFAULT_OP = 0  // 0 = original OpenMesh implementation
  };


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

  void set_lock(bool lock) { lock_boundary_edges = lock; }

  void set_min_mod(int mod_) {
    switch(mod_){
      case 0: mod = DEFAULT_OP; break;
      case 1: mod = POINTS;     break;  
      case 2: mod = LINE;       break;  
      case 3: mod = SPACE;      break;
      default: mod = DEFAULT_OP;
    }
   }

  void set_opts(std::string opts) {

    //parse 1st option (lock boundary edges)
    size_t pos = opts.find(",");
    bool lock;
    std::string first = opts.substr(0, pos);            
    if (first == "true" or first == "false") {
      std::istringstream(first) >> std::boolalpha >> lock;
      set_lock(lock);
    }
    else if (first.size() == 0) {}
    else std::cerr << "Invalid first option - either \"true\" or \"false\" required, default value (false) was set." << std::endl;

    if (pos != std::string::npos) {
      opts.erase(0, pos+1);
      pos = opts.find(",");
      std::string first = opts.substr(0, pos);
      int number = std::stoi(first);
      if (number == 0 or number == 1 or number == 2 or number == 3) set_min_mod(number);
      else {
        std::cout<<"Invalid second option - permitted values: 0, 1, 2, 3\n"<<
                  "0 = OpenMesh implementation; 1 = start/end/mid points only;\
                  2 = line determined by v0 and v1; 3 = original GH (anywhere in 3D space)" << std::endl;
        std::cout<<"Going with original OpenMesh implementation"<<std::endl;
      }  
      
    }

    if (pos != std::string::npos) {
      opts.erase(0, pos+1);
      set_max_err(std::stoi(opts));
    }

  }


private:

   // ------Parameters of the decimating module---------------------------------


  /** Parameter for locking all boundary and "semi-boundary" edge (edges with 
   * only one boundary vertex) to preserve the mesh boundary. If we don't
   * collapse these edges, the boundary will stay the same.
   */
  bool lock_boundary_edges = false;

  // maximum quadric error
  double max_err_;

  // Decimation mode parameter 
  decimation_mode mod = DEFAULT_OP;

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

  /** Ideal vertex position after collapse
   * Solution to equation Av=b
   * solved using: v = A^(-1) * b
   */
  HPropHandleT<DefaultTraits::Point> ideal_vertex_coords;

  // this vertex property stores a quadric for each vertex
  VPropHandleT< Geometry::QuadricT<double> >  quadrics_;
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

