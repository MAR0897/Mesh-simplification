//=============================================================================
//
//  CLASS ModSpectralT
//
//=============================================================================
#ifndef OSG_MODSPECTRAL_HH
#define OSG_MODSPECTRAL_HH
//== INCLUDES =================================================================
#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <set>

#include <OpenMesh/Eigen/Dense>
#include <OpenMesh/Eigen/Sparse>
#include <OpenMesh/Spectra/include/Spectra/SymEigsShiftSolver.h>
#include <OpenMesh/Spectra/include/Spectra/MatOp/SparseSymShiftSolve.h>


//== NAMESPACE ================================================================
namespace OpenMesh  {
namespace Decimater {
//== CLASS DEFINITION =========================================================


/** \brief Mesh decimation module computing collapse priority based on .
 *
 *  
 */
template <class MeshT>
class ModSpectralT : public ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModSpectralT, MeshT, Spectral );

  //Matrix is an array of rows (which is Vec3d)
  using Matrix3d = std::array<Vec3d, 3>;

public:

   explicit ModSpectralT( MeshT &_mesh )
    : Base(_mesh, false)
  {
    // Add needed mesh properties for Lind-Turk decimation
    Base::mesh().add_property(SPprops);
    Base::mesh().add_property(idx);
    Base::mesh().add_property(local_idx);
    Base::mesh().add_property(area);
    Base::mesh().add_property(cotangents);
  }
  
  virtual ~ModSpectralT()
  {
    Base::mesh().remove_property(SPprops);
    Base::mesh().remove_property(idx);
    Base::mesh().remove_property(local_idx);
    Base::mesh().remove_property(area);
    Base::mesh().remove_property(cotangents);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation, possibly lock boundary edges if option is set
  virtual void initialize(void) override;

  // Compute error and remaining vertex position for a halfedge
  virtual float collapse_priority(const CollapseInfo& _ci) override;
	
  //tu posun v1 na pozici p1 pomoci set_point
  virtual void preprocess_collapse(const CollapseInfo& _ci) override;

  //update signals and previous error
  virtual void postprocess_collapse(const CollapseInfo& _ci) override;
 
  void set_eigenvec_n(const double& _eigenvec_n) { eigenvec_n = _eigenvec_n; }

  void set_lock(bool lock) { lock_boundary_edges = lock; }

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
      if (number > 0) set_eigenvec_n(number);  
    }
  }

  double calc_cotangent(const HalfedgeHandle& he, const VertexHandle& vh1, const VertexHandle& vh2);
  //inner function, p3 is the vertex where the angle is being computed
  double calc_cotangent_from_points(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3);

private:

  //command line args
  bool lock_boundary_edges = false;
  int eigenvec_n = 100;

	//matrices
	Eigen::MatrixXd F;              //signals
	Eigen::MatrixXd Z;              //signals_L
	
  //cost
  Eigen::VectorXd norms;  //previous error
	
  //Properties
  struct Props {
    bool is_locked;           //optional lock for boundary and 'semi-boundary' vertices to preserve mesh boundary
    bool error_calculated;
    double alpha;             //collapse polynomial minimum
    Vec3d res_vertex_coords;  //ideal resulting vertex for collapsed edge     
    std::vector<std::pair<VertexHandle, double>> cost_diff; //cost difference for each 1-ring vertex of the KEEP vertex     
    //std::set<FaceHandle> recalc_faces;  //faces to recalculate area and cotangents for
  };

  HPropHandleT<Props>  SPprops;
  VPropHandleT<size_t> idx;
  VPropHandleT<size_t> local_idx;
  FPropHandleT<double> area;
  FPropHandleT<std::unordered_map<VertexHandle, double>> cotangents;
};

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_DECIMATER_MODSPECTRAL_CC)
#define OSG_MODSPECTRAL_TEMPLATES
#include "ModSpectralT_impl.hh"
#endif
//=============================================================================
#endif // OSG_MODLINDTURK_HH defined
//=============================================================================
