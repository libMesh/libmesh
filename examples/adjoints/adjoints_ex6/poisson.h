#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class PoissonSystem : public FEMSystem
{
public:
  PoissonSystem(EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _fe_family("LAGRANGE"), _fe_order(2),
      _analytic_jacobians(true) { qoi.resize(1); computed_QoI[0] = 0.0; }

  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // Postprocessing function which we are going to override for this application

  virtual void postprocess(void);

  Number & get_QoI_value(std::string type, unsigned int QoI_index)
  {
    if (type == "exact")
      {
        return exact_QoI[QoI_index];
      }
    else
      {
        return computed_QoI[QoI_index];
      }
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

  // Overloading the postprocess function

  virtual void element_postprocess(DiffContext & context);

  // Forcing function coefficient
  Real alpha;

  // Indices for each variable
  unsigned int T_var;

  // Variables to hold the computed QoIs

  Number computed_QoI[1];

  // Variables to read in the exact QoIs from poisson.in

  Number exact_QoI[1];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
