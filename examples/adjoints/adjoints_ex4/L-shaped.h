#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"

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
                const unsigned int number_in) :
    FEMSystem(es, name_in, number_in),
    _fe_family("LAGRANGE"),
    _fe_order(1),
    _analytic_jacobians(true)
  { qoi.resize(2); }

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // Postprocessing function which we are going to override for this application

  virtual void postprocess();

  Number & get_QoI_value(std::string type, unsigned int QoI_index)
  {
    if (type == "exact")
      return exact_QoI[QoI_index];
    else
      return computed_QoI[QoI_index];
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

  // Overloading the postprocess function

  virtual void element_postprocess(DiffContext & context);

  virtual void side_postprocess(DiffContext & context);

  // Overloading the qoi function on elements

  virtual void element_qoi_derivative (DiffContext & context,
                                       const QoISet & qois);

  // Overloading the qoi function on sides

  virtual void side_qoi_derivative (DiffContext & context,
                                    const QoISet & qois);

  Number exact_solution (const Point &);

  // Variables to hold the computed QoIs

  Number computed_QoI[2];

  // Variables to read in the exact QoIs from l-shaped.in

  Number exact_QoI[2];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
