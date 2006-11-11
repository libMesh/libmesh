#include "diff_solver.h"
#include "diff_system.h"
#include "steady_solver.h"



SteadySolver::~SteadySolver ()
{
}
  


bool SteadySolver::element_residual(bool request_jacobian)
{
  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // The user shouldn't compute a jacobian unless requested
  assert(request_jacobian || !jacobian_computed);

  bool jacobian_computed2 =
    _system.element_constraint(jacobian_computed);

  // User code should either always return an analytic jacobian
  // or never return one
  assert(jacobian_computed == jacobian_computed2);

  return jacobian_computed && jacobian_computed2;
}



bool SteadySolver::side_residual(bool request_jacobian)
{
  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  // The user shouldn't compute a jacobian unless requested
  assert (request_jacobian || !jacobian_computed);

  bool jacobian_computed2 =
    _system.side_constraint(jacobian_computed);

  // User code should either always return an analytic jacobian
  // or never return one
  assert (jacobian_computed == jacobian_computed2);
  
  return jacobian_computed && jacobian_computed2;
}

