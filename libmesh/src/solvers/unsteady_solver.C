// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include "diff_solver.h"
#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "unsteady_solver.h"



UnsteadySolver::UnsteadySolver (sys_type& s)
  : TimeSolver(s),
    old_local_nonlinear_solution (NumericVector<Number>::build()),
    first_solve                  (true)
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

#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_nonlinear_solution->init (_system.n_dofs(), _system.n_local_dofs(),
                                      _system.get_dof_map().get_send_list(), false,
                                      GHOSTED);
#else
  old_local_nonlinear_solution->init (_system.n_dofs(), false, SERIAL);
#endif

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
	      libMesh::out << "Newton backtracking failed.  Trying with smaller timestep, dt="
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
		    libMesh::out << "Reduced dt solve succeeded." << std::endl;
		  return;
		}
	    }

	  // If we made it here, we still couldn't converge the solve after
	  // reducing deltat
	  libMesh::out << "DiffSolver::solve() did not succeed after "
		        << reduce_deltat_on_diffsolver_failure
		        << " attempts." << std::endl;
	  libmesh_convergence_failure();
	  
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
  libmesh_assert (global_dof_number < _system.get_dof_map().n_dofs());
  libmesh_assert (global_dof_number < old_local_nonlinear_solution->size());

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
