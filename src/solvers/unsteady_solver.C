// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/unsteady_solver.h"

namespace libMesh
{



UnsteadySolver::UnsteadySolver (sys_type & s)
  : TimeSolver(s),
    old_local_nonlinear_solution (NumericVector<Number>::build(s.comm())),
    first_solve                  (true),
    first_adjoint_step (true)
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



void UnsteadySolver::init_data()
{
  TimeSolver::init_data();

#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_nonlinear_solution->init (_system.n_dofs(), _system.n_local_dofs(),
                                      _system.get_dof_map().get_send_list(), false,
                                      GHOSTED);
#else
  old_local_nonlinear_solution->init (_system.n_dofs(), false, SERIAL);
#endif
}



void UnsteadySolver::reinit ()
{
  TimeSolver::reinit();

#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_nonlinear_solution->init (_system.n_dofs(), _system.n_local_dofs(),
                                      _system.get_dof_map().get_send_list(), false,
                                      GHOSTED);
#else
  old_local_nonlinear_solution->init (_system.n_dofs(), false, SERIAL);
#endif

  // localize the old solution
  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");

  old_nonlinear_soln.localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}



void UnsteadySolver::solve ()
{
  if (first_solve)
    {
      advance_timestep();
      first_solve = false;
    }

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
              bool backtracking_still_failed =
                solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;

              bool backtracking_max_iterations =
                solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;

              if (!backtracking_still_failed && !backtracking_max_iterations)
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
  // The first access of advance_timestep happens via solve, not user code
  // It is used here to store any initial conditions data
  if (!first_solve)
    {
      // We call advance_timestep in user code after solve, so any solutions
      // we will be storing will be for the next time instance
      _system.time += _system.deltat;
    }
    else
    {
      // We are here because of a call to advance_timestep that happens
      // via solve, the very first solve. All we are doing here is storing
      // the initial condition. The actual solution computed via this solve
      // will be stored when we call advance_timestep in the user's timestep loop
      first_solve = false;
    }

  // If the user has attached a memory or file solution history object
  // to the solver, this will store the current solution indexed with
  // the current time
  solution_history->store(false);

  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> & nonlinear_solution =
    *(_system.solution);

  old_nonlinear_soln = nonlinear_solution;

  old_nonlinear_soln.localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}



void UnsteadySolver::adjoint_advance_timestep ()
{
  // All calls to adjoint_advance_timestep are made in the user's
  // code. This first call is made immediately after the adjoint initial conditions
  // are set. This is in the user code outside the adjoint time loop.
  if (!first_adjoint_step)
    {
      // The adjoint system has been solved. We need to store the adjoint solution and
      // load the primal solutions for the next time instance (t - delta_ti).
      _system.time -= _system.deltat;
    }
  else
    {
      // The first adjoint step simply saves the given adjoint initial condition
      // So there is a store, but no solve, no actual timestep, so no need to change system time
      first_adjoint_step = false;
    }

  // Retrieve the primal solution vectors at this new (or for
  // first_adjoint_step, initial) time instance. These provide the
  // data to solve the adjoint problem for the next time instance.
  solution_history->retrieve(true);

  // Call the store function to store the adjoint we have computed (or
  // for first_adjoint_step, the adjoint initial condition) in this
  // time step for the time instance.
  solution_history->store(true);

  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

void UnsteadySolver::retrieve_timestep()
{
  // Retrieve all the stored vectors at the current time
  solution_history->retrieve(false);

  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}


Number UnsteadySolver::old_nonlinear_solution(const dof_id_type global_dof_number)
  const
{
  libmesh_assert_less (global_dof_number, _system.get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, old_local_nonlinear_solution->size());

  return (*old_local_nonlinear_solution)(global_dof_number);
}



Real UnsteadySolver::du(const SystemNorm & norm) const
{

  std::unique_ptr<NumericVector<Number>> solution_copy =
    _system.solution->clone();

  solution_copy->add(-1., _system.get_vector("_old_nonlinear_solution"));

  solution_copy->close();

  return _system.calculate_norm(*solution_copy, norm);
}

} // namespace libMesh
