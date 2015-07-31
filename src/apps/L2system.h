#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class L2System : public libMesh::FEMSystem
{
public:
  // Constructor
  L2System(libMesh::EquationSystems& es,
           const std::string& name,
           const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    _fe_family("LAGRANGE"), _fe_order(1) {}

  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }
  libMesh::AutoPtr<libMesh::FEMFunctionBase<libMesh::Number> > goal_func;

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (libMesh::DiffContext &context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext &context);

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;
};
