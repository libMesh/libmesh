
#include "diff_solver.h"
#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "unsteady_solver.h"



UnsteadySolver::UnsteadySolver (sys_type& s)
  : TimeSolver(s),
    first_solve                  (true),
    old_local_nonlinear_solution (NumericVector<Number>::build())
{
}



UnsteadySolver::~UnsteadySolver ()
{
}



void UnsteadySolver::init ()
{
  TimeSolver::init();

  _system.add_vector("_old_nonlinear_solution");
}



void UnsteadySolver::solve ()
{
  if (first_solve)
    {
      advance_timestep();
      first_solve = false;
    }

  old_local_nonlinear_solution->init (_system.n_dofs());

  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());

  unsigned int solve_result = _diff_solver->solve();

  // If we requested the UnsteadySolver to attempt reducing dt after a
  // failed DiffSolver solve, check the results of the solve now.
  if (reduce_deltat_on_diffsolver_failure)
    {
      bool backtracking_failed =
	solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;

      bool max_iterations =
	solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
	
      if (backtracking_failed || max_iterations)
	{
	  // Cut timestep in half
	  for (unsigned int nr=0; nr<reduce_deltat_on_diffsolver_failure; ++nr)
	    {
	      _system.deltat *= 0.5;
	      std::cout << "Newton backtracking failed.  Trying with smaller timestep, dt="
			<< _system.deltat << std::endl;

	      solve_result = _diff_solver->solve();

	      // Check solve results with reduced timestep
	      bool backtracking_failed =
		solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;
	      
	      bool max_iterations =
		solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;

	      if (!backtracking_failed && !max_iterations)
		{
		  if (!quiet)
		    std::cout << "Reduced dt solve succeeded." << std::endl;
		  return;
		}
	    }

	  // If we made it here, we still couldn't converge the solve after
	  // reducing deltat
	  std::cout << "DiffSolver::solve() did not succeed after "
		    << reduce_deltat_on_diffsolver_failure
		    << " attempts." << std::endl;
	  error();
	  
  	} // end if (backtracking_failed || max_iterations)
    } // end if (reduce_deltat_on_diffsolver_failure)
}



void UnsteadySolver::advance_timestep ()
{
  NumericVector<Number> &old_nonlinear_solution =
  _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> &nonlinear_solution =
    *(_system.solution);

  old_nonlinear_solution = nonlinear_solution;

  if (!first_solve)
    _system.time += _system.deltat;
}



Number UnsteadySolver::old_nonlinear_solution(const unsigned int global_dof_number)
const
{
  assert (global_dof_number < _system.get_dof_map().n_dofs());
  assert (global_dof_number < old_local_nonlinear_solution->size());

  return (*old_local_nonlinear_solution)(global_dof_number);
}



Real UnsteadySolver::du(const SystemNorm &norm) const
{

  AutoPtr<NumericVector<Number> > solution_copy =
    _system.solution->clone();

  solution_copy->add(-1., _system.get_vector("_old_nonlinear_solution"));

  solution_copy->close();

  return _system.calculate_norm(*solution_copy, norm);
}
