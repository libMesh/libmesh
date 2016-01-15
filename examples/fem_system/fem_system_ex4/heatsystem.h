#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class HeatSystem : public libMesh::FEMSystem
{
public:
  // Constructor
  HeatSystem(libMesh::EquationSystems & es,
             const std::string & name,
             const unsigned int number) :
    libMesh::FEMSystem(es, name, number),
    _fe_family("LAGRANGE"), _fe_order(1)
  {
    // Get the conductivity ratios right for both 2D and 3D
    // benchmarks
    _k[1] = 1/libMesh::pi/std::sqrt(libMesh::Real(3));
    _k[2] = 1;
    _k[3] = 2*libMesh::pi*std::sqrt(libMesh::Real(3));
  }

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context);

  // The conductivity for the various dimensional elements, indexed by
  // dim (with _k[0] unused) for simplicity
  libMesh::Real _k[4];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // The variable index (yes, this will be 0...)
  unsigned int T_var;
};
