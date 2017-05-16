#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"

#include <map>

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class L2System : public libMesh::FEMSystem
{
public:
  // Constructor
  L2System(libMesh::EquationSystems & es,
           const std::string & name,
           const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    input_system(NULL),
    _fe_family("LAGRANGE"),
    _fe_order(1) {}

  // Destructor; deletes extra context objects
  ~L2System();

  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }

  // We want to be able to project functions based on *other* systems'
  // values.  For that we need not only a FEMFunction but also a
  // reference to the system where it applies and a separate context
  // object (or multiple separate context objects, in the threaded
  // case) for that system.
  libMesh::AutoPtr<libMesh::FEMFunctionBase<libMesh::Number> > goal_func;

  libMesh::System * input_system;

  std::map<libMesh::FEMContext *, libMesh::FEMContext *>
    input_contexts;

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context);

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;
};
