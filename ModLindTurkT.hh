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

  //Matrix is an array of rows (which is Vec3d)
  using Matrix3d = std::array<Vec3d, 3>;

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
	
  //ZATIM NIC NEZADAVEJ, ALE KDYZ TAK SE K TOMUTO BODU VRAT: tady do CollapseInfo zadej pozici nejlepsiho umisteni remaining vertexu = p1
  //PRAVDEPODOBNE TO COLLAPSE INFO JE JEN INFO K VYPOCTU ERRORU, NIC SE PODLE TOHO NEKOLABUJE, TAKZE NEMUSIS NIC VPISOVAT
  //tu posun v1 na pozici p1 pomoci set_point
  virtual void preprocess_collapse(const CollapseInfo& _ci) override;
    
  //a nebo tu posun v1 na pozici p1 pomoci set_point, melo by to byt jedno, ta kvadricka simplifikace to stejne vubec neposunuje
  //virtual void postprocess_collapse(const CollapseInfo& _ci) override;

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
    opts.erase(0, pos+1);

  
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
  
  void set_lambda(double _lambda) { lambda = _lambda; }

  void set_lock(double _lock) { lock_boundary_edges = _lock; }

  void set_alpha(double _alpha = 0.01745329251)
  {
    alpha = _alpha;
    SINALPHA2 = std::pow(std::sin(alpha), 2);
    COSALPHA2 = std::pow(std::cos(alpha), 2);
  }

  //Checks if constraint is alpha compatible, so we can add it to the system using 'add_constraint'
  bool is_alpha_compatible(const HalfedgeHandle& heh, const Vec3d& constraint);

  //Adds constraint to current halfedge properties
  void add_constraint(const HalfedgeHandle& heh, const Vec3d& constraint, double& right_side);

  //Calculates remaining constraints if there are less than 3
  void calc_remaining_constraints(const HalfedgeHandle& heh, Matrix3d& Hessian, Vec3d& c);

  Matrix3d inverse3x3(const Matrix3d& m);

private:

  bool A = true;
  //std::vector<HalfedgeHandle> error;
  bool lock_boundary_edges = false;
  double  lambda = 0.5,
          alpha,
          SINALPHA2,
          COSALPHA2;

  struct Props {
    bool is_locked;           //optional lock for boundary and 'semi-boundary' vertices to preserve mesh boundary
    size_t n;                 //number of valid constraints 
    //double cost;              //the cost (error) of collapsing the edge
    Vec3d res_vertex_coords;  //ideal resulting vertex for collapsed edge
    Vec3d b_side;             //b side vector of the system of constraints 
    Matrix3d constraints;     //constraining the resulting vertex to one point               
    //int ncalc;
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
