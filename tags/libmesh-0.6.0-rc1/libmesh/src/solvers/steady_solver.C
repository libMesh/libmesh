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

  jacobian_computed =
    _system.element_constraint(request_jacobian) && jacobian_computed;

  return jacobian_computed;
}



bool SteadySolver::side_residual(bool request_jacobian)
{
  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  jacobian_computed =
    _system.side_constraint(request_jacobian) && jacobian_computed;
  
  return jacobian_computed;
}

