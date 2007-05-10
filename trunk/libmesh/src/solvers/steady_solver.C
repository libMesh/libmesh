#include "diff_solver.h"
#include "diff_system.h"
#include "steady_solver.h"



SteadySolver::~SteadySolver ()
{
}
  


bool SteadySolver::element_residual(bool request_jacobian)
{
  // If a fixed solution is requested, it will just be the current
  // solution
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = _system.elem_solution;
      _system.fixed_solution_derivative = 1.0;
    }

  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // The user shouldn't compute a jacobian unless requested
  assert(request_jacobian || !jacobian_computed);

  bool jacobian_computed2 =
    _system.element_constraint(jacobian_computed);

  // The user shouldn't compute a jacobian unless requested
  assert (jacobian_computed || !jacobian_computed2);

  return jacobian_computed && jacobian_computed2;
}



bool SteadySolver::side_residual(bool request_jacobian)
{
  // If a fixed solution is requested, it will just be the current
  // solution
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = _system.elem_solution;
      _system.fixed_solution_derivative = 1.0;
    }

  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  // The user shouldn't compute a jacobian unless requested
  assert (request_jacobian || !jacobian_computed);

  bool jacobian_computed2 =
    _system.side_constraint(jacobian_computed);

  // The user shouldn't compute a jacobian unless requested
  assert (jacobian_computed || !jacobian_computed2);
  
  return jacobian_computed && jacobian_computed2;
}

