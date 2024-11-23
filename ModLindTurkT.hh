//=============================================================================
//
//  CLASS ModLindTurkT
//
//=============================================================================
#ifndef OSG_MODLINDTURK_HH
#define OSG_MODLINDTURK_HH
//== INCLUDES =================================================================
#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Eigen/Dense>
//== NAMESPACE ================================================================
namespace OpenMesh  {
namespace Decimater {
//== CLASS DEFINITION =========================================================


/** \brief Mesh decimation module computing collapse priority based on .
 *
 *  
 */
template <class MeshT>
class ModLindTurkT : public ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModLindTurkT, MeshT, LindTurk );

public:

   explicit ModLindTurkT( MeshT &_mesh )
    : Base(_mesh, false)
  {
    // Add needed mesh properties for Lind-Turk decimation
    Base::mesh().add_property(LTprops);
  }
  
  virtual ~ModLindTurkT()
  {
    Base::mesh().remove_property(LTprops);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation, possibly lock boundary edges if option is set
  virtual void initialize(void) override;

  // Compute error and remaining vertex position for a halfedge
  virtual float collapse_priority(const CollapseInfo& _ci) override;
	
  //TO COLLAPSE INFO JE JEN INFO K VYPOCTU ERRORU, NIC SE PODLE TOHO NEKOLABUJE, TAKZE NETREBA NIC VPISOVAT
  //tu posun v1 na pozici p1 pomoci set_point
  virtual void preprocess_collapse(const CollapseInfo& _ci) override {
    //move remaining vertex to ideal calculated position
    Base::mesh().set_point(_ci.v1, Base::mesh().property(LTprops, _ci.v0v1).res_vertex_coords);
  }

  void set_opts(std::string opts)
  {
     //parse 1st option (lock boundary edges)
    size_t pos = opts.find(",");
    bool lock;
    std::string first = opts.substr(0, pos);            
    if (first == "true" or first == "false") {
      std::istringstream(first) >> std::boolalpha >> lock;
      set_lock(lock);
    }
    else if (first.size() == 0) {}
    else std::cerr << "Invalid first LT option - either \"true\" or \"false\" required, default value (false) was set." << std::endl;

    if (pos != std::string::npos) {
      opts.erase(0, pos+1);
      std::cout<<"String: "<<opts<<std::endl;
      //parse 2nd option (lambda - final error calculation weight)
      if (pos != std::string::npos) {
        size_t pos = opts.find(",");
        auto second = opts.substr(0, pos);
        if (second.size() != 0) set_lambda(std::stod(second));
        opts.erase(0, pos+1); 
      }
      //parse 3rd option (alpha - angle to which planes are taken as coplanar)
      if (!opts.empty()) {
        set_alpha(std::stod(opts));
        opts.erase(0, pos+1); 
      }
    }
    else set_alpha();
  }
  
  void set_lambda(double _lambda) { lambda = _lambda; }

  void set_lock(double _lock) { lock_boundary_edges = _lock; }

  void set_alpha(double _alpha = 0.01745329251)
  {
    alpha = _alpha;
    SINALPHA2 = std::pow(std::sin(alpha), 2);
    COSALPHA2 = std::pow(std::cos(alpha), 2);
  }

  //Checks if constraint is alpha compatible, so we can add it to the system using 'add_constraint'
  bool is_alpha_compatible(const HalfedgeHandle& heh, const Eigen::Vector3d& constraint);

  //Adds constraint to current halfedge properties
  void add_constraint(const HalfedgeHandle& heh, const Eigen::Vector3d& constraint, const double& right_side);

  //Calculates remaining constraints if there are less than 3
  void calc_remaining_constraints(const HalfedgeHandle& heh, const Eigen::Matrix3d& Hessian, const Eigen::Vector3d& c);

private:

  bool lock_boundary_edges = false;
  double  lambda = 0.5,
          alpha,
          SINALPHA2,
          COSALPHA2;

  struct Props {
    //checks if the cost error of a halfedge has been computed,
    // so we dont need to compute it for the other halfedge,
    // as this algorithm only work with edges 
    bool error_calculated;    
    bool is_locked;           //optional lock for boundary and 'semi-boundary' vertices to preserve mesh boundary
    size_t n;                 //number of valid constraints 
    DefaultTraits::Point res_vertex_coords;  //ideal resulting vertex for collapsed edge
    Eigen::Vector3d b_side;             //b side vector of the system of constraints 
    Eigen::Matrix3d constraints;     //constraining the resulting vertex to one point               
  };

  HPropHandleT<Props>  LTprops;
};

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_DECIMATER_MODLINDTURK_CC)
#define OSG_MODLINDTURK_TEMPLATES
#include "ModLindTurkT_impl.hh"
#endif
//=============================================================================
#endif // OSG_MODLINDTURK_HH defined
//=============================================================================
