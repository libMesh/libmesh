// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/numeric_vector.h"

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
  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> & nonlinear_solution =
    *(_system.solution);
  //    _system.get_vector("_nonlinear_solution");

  old_nonlinear_soln = nonlinear_solution;

  if (!first_solve)
    _system.time += last_deltat;
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
