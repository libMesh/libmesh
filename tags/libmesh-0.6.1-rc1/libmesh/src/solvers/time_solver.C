
#include "diff_solver.h"
#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s)
  : quiet(true),
    reduce_deltat_on_diffsolver_failure(0),
    _diff_solver                 (NULL),
    _system                      (s),
    first_solve                  (true),
    old_local_nonlinear_solution (NumericVector<Number>::build())
{
}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::reinit ()
{
  _diff_solver->reinit();
}



void TimeSolver::init ()
{
  // If the user hasn't given us a solver to use,
  // just build a default solver
  if (_diff_solver.get() == NULL)
    _diff_solver = DiffSolver::build(_system);
  _diff_solver->init();

  _system.add_vector("_old_nonlinear_solution");
}



void TimeSolver::solve ()
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

  // If we requested the TimeSolver to attempt reducing dt after a
  // failed DiffSolver solve, check the results of the solve now.
  if (reduce_deltat_on_diffsolver_failure)
    {
      bool backtracking_failed =
	solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;

      bool max_iterations =
	solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
	
      if (backtracking_failed || max_iterations)
	{
	  // If the solve failed for one of the reasons above and this is
	  // a SteadySolver, then we can't try to make it succeed by
	  // decreasing deltat.  We can test to see if we are a
	  // SteadySolver by checking our error_order().
	  if  (this->error_order() == 0.0)
	    {
	      std::cout << "DiffSolver::solve() failed, but we cannot reduce deltat "
			<< " because this is a SteadySolver!" << std::endl;
	      error();
	    }

	  // Otherwise, cut timestep in half
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



void TimeSolver::advance_timestep ()
{
  NumericVector<Number> &old_nonlinear_solution =
  _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> &nonlinear_solution =
    *(_system.solution);
//    _system.get_vector("_nonlinear_solution");

  old_nonlinear_solution = nonlinear_solution;

  if (!first_solve)
    _system.time += _system.deltat;
}



Number TimeSolver::old_nonlinear_solution(const unsigned int global_dof_number)
const
{
  assert (global_dof_number < _system.get_dof_map().n_dofs());
  assert (global_dof_number < old_local_nonlinear_solution->size());

  return (*old_local_nonlinear_solution)(global_dof_number);
}



AutoPtr<DiffSolver> & TimeSolver::diff_solver()
{
  return _diff_solver;
}



Real TimeSolver::du(unsigned char norm_type)
{

  AutoPtr<NumericVector<Number> > solution_copy =
    _system.solution->clone();

  solution_copy->add(-1., _system.get_vector("_old_nonlinear_solution"));

  solution_copy->close();

  if (norm_type==0)
    return solution_copy->l2_norm();

  else if (norm_type==1)
    return solution_copy->l1_norm();

  else
    {
      std::cout << "Unrecognized norm "
		<< norm_type
		<< " requested.";
      error();
    }

  return 0.;
    
}
