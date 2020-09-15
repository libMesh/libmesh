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

#include "libmesh/adaptive_time_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/no_solution_history.h"

namespace libMesh
{



AdaptiveTimeSolver::AdaptiveTimeSolver (sys_type & s)
  : FirstOrderUnsteadySolver(s),
    core_time_solver(),
    target_tolerance(1.e-3),
    upper_tolerance(0.0),
    max_deltat(0.),
    min_deltat(0.),
    max_growth(0.),
    completed_timestep_size(1.),
    global_tolerance(true)
{
  // the child class must populate core_time_solver
  // with whatever actual time solver is to be used

  // As an UnsteadySolver, we have an old_local_nonlinear_solution, but we're
  // going to drop it and use our core time solver's instead.
  old_local_nonlinear_solution.reset();
}



AdaptiveTimeSolver::~AdaptiveTimeSolver ()
{
  // As an UnsteadySolver, we have an old_local_nonlinear_solution, but it
  // is being managed by our core_time_solver.  Make sure we don't delete
  // it out from under them, in case the user wants to keep using the core
  // solver after they're done with us.
  old_local_nonlinear_solution.release();
}



void AdaptiveTimeSolver::init()
{
  libmesh_assert(core_time_solver.get());

  // We override this because our core_time_solver is the one that
  // needs to handle new vectors, diff_solver->init(), etc
  core_time_solver->init();

  // Set the core_time_solver's solution history object to be the same one as
  // that for the outer adaptive time solver
  core_time_solver->set_solution_history((this->get_solution_history()));

  // Now that we have set the SolutionHistory object for the coretimesolver,
  // we set the SolutionHistory type for the timesolver to be NoSolutionHistory
  // All storage and retrieval will be handled by the coretimesolver directly.
  NoSolutionHistory outersolver_solution_history;
  this->set_solution_history(outersolver_solution_history);

  // As an UnsteadySolver, we have an old_local_nonlinear_solution, but it
  // isn't pointing to the right place - fix it
  //
  // This leaves us with two std::unique_ptrs holding the same pointer - dangerous
  // for future use.  Replace with shared_ptr?
  old_local_nonlinear_solution =
    std::unique_ptr<NumericVector<Number>>(core_time_solver->old_local_nonlinear_solution.get());
}



void AdaptiveTimeSolver::reinit()
{
  libmesh_assert(core_time_solver.get());

  // We override this because our core_time_solver is the one that
  // needs to handle new vectors, diff_solver->reinit(), etc
  core_time_solver->reinit();
}


void AdaptiveTimeSolver::advance_timestep ()
{
  // The first access of advance_timestep happens via solve, not user code
  // It is used here to store any initial conditions data
  if (!first_solve)
    {
      _system.time += this->completed_timestep_size;
    }
  else
    {
    // We are here because of a call to advance_timestep that happens
    // via solve, the very first solve. All we are doing here is storing
    // the initial condition. The actual solution computed via this solve
    // will be stored when we call advance_timestep in the user's timestep loop
    first_solve = false;
    }

  // For the adaptive time solver, all SH operations
  // are handled by the core_time_solver's SH object
  // Sub solution storage is handled internally by the core time solver,
  // but the 'full step' solution is stored here to maintain consistency
  // with the fixed timestep scheme.
  core_time_solver->get_solution_history().store(false, _system.time);

  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> & nonlinear_solution =
    *(_system.solution);

  old_nonlinear_soln = nonlinear_solution;

  old_nonlinear_soln.localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

void AdaptiveTimeSolver::adjoint_advance_timestep ()
{
  // All calls to adjoint_advance_timestep are made in the user's
  // code. This first call is made immediately after the adjoint initial conditions
  // are set. This is in the user code outside the adjoint time loop.
  if (!first_adjoint_step)
    {
      // The adjoint system has been solved. We need to store the adjoint solution and
      // load the primal solutions for the next time instance (t - delta_ti).
      _system.time -= this->completed_timestep_size;
    }
  else
    {
      // The first adjoint step simply saves the given adjoint initial condition
      // So there is a store, but no solve, no actual timestep, so no need to change system time
      first_adjoint_step = false;
    }

  // For the adaptive time solver, all SH operations
  // are handled by the core_time_solver's SH object
  // Retrieve the primal solution for the next adjoint calculation,
  // by using the core time solver's solution history object.
  core_time_solver->get_solution_history().retrieve(true, _system.time);

  // Store the computed full step adjoint solution for future use (sub steps are handled internally by the core time solver)
  core_time_solver->get_solution_history().store(true, _system.time);

  // We also need to tell the core time solver that the adjoint initial conditions have been set
  core_time_solver->set_first_adjoint_step(false);

  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
    _system.get_dof_map().get_send_list());
}

void AdaptiveTimeSolver::retrieve_timestep ()
{
  // Ask the core time solver to retrieve all the stored vectors
  // at the current time
  core_time_solver->retrieve_timestep();
}

Real AdaptiveTimeSolver::error_order () const
{
  libmesh_assert(core_time_solver.get());

  return core_time_solver->error_order();
}



bool AdaptiveTimeSolver::element_residual (bool request_jacobian,
                                           DiffContext & context)
{
  libmesh_assert(core_time_solver.get());

  return core_time_solver->element_residual(request_jacobian, context);
}



bool AdaptiveTimeSolver::side_residual (bool request_jacobian,
                                        DiffContext & context)
{
  libmesh_assert(core_time_solver.get());

  return core_time_solver->side_residual(request_jacobian, context);
}



bool AdaptiveTimeSolver::nonlocal_residual (bool request_jacobian,
                                            DiffContext & context)
{
  libmesh_assert(core_time_solver.get());

  return core_time_solver->nonlocal_residual(request_jacobian, context);
}



std::unique_ptr<DiffSolver> & AdaptiveTimeSolver::diff_solver()
{
  return core_time_solver->diff_solver();
}



std::unique_ptr<LinearSolver<Number>> & AdaptiveTimeSolver::linear_solver()
{
  return core_time_solver->linear_solver();
}



Real AdaptiveTimeSolver::calculate_norm(System & s,
                                        NumericVector<Number> & v)
{
  return s.calculate_norm(v, component_norm);
}

} // namespace libMesh
