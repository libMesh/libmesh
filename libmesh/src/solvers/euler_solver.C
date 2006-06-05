
#include "euler_solver.h"
#include "diff_system.h"



EulerSolver::EulerSolver (sys_type& s)
 : TimeSolver(s), theta(1.)
{
}



EulerSolver::~EulerSolver ()
{
}



void EulerSolver::init ()
{
  Parent::init();

  _system.add_vector("_old_nonlinear_solution");
}



bool EulerSolver::element_residual (bool request_jacobian)
{
  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  jacobian_computed = _system.element_constraint(request_jacobian) &&
    jacobian_computed;

  return jacobian_computed;
}



bool EulerSolver::side_residual (bool request_jacobian)
{
  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  jacobian_computed = _system.side_constraint(request_jacobian) &&
    jacobian_computed;

  return jacobian_computed;
}
