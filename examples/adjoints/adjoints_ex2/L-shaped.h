#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"
#include "libmesh/parameter_pointer.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class LaplaceSystem : public FEMSystem
{
public:
  // Constructor
  LaplaceSystem(EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _fe_family("LAGRANGE"),
      _fe_order(1),
      _analytic_jacobians(true)
  {}

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  Number & get_parameter_value(unsigned int parameter_index)
  {
    return parameters[parameter_index];
  }

  ParameterVector & get_parameter_vector()
  {
    typedef ParameterPointer<Number> PP;
    parameter_vector.clear();
    for (std::size_t i = 0; i != parameters.size(); ++i)
      parameter_vector.push_back(UniquePtr<ParameterAccessor<Number> >(new PP(&parameters[i])));

    return parameter_vector;
  }

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Constraint parts
  virtual bool side_constraint (bool request_jacobian,
                                DiffContext & context);

  Number exact_solution (const Point &);

  // Parameters associated with the system
  std::vector<Number> parameters;

  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
