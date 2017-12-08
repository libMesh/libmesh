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


#include "libmesh/libmesh_common.h"
#include "libmesh/diff_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/implicit_system.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

DiffSolver::DiffSolver (sys_type & s) :
  ParallelObject(s),
  max_linear_iterations(1000),
  max_nonlinear_iterations(100),
  quiet(true),
  verbose(false),
  continue_after_max_iterations(true),
  continue_after_backtrack_failure(false),
  absolute_residual_tolerance(0.),
  relative_residual_tolerance(0.),
  absolute_step_tolerance(0.),
  relative_step_tolerance(0.),
  initial_linear_tolerance(1e-12),
  minimum_linear_tolerance(TOLERANCE*TOLERANCE),
  max_solution_norm(0.),
  max_residual_norm(0.),
  _outer_iterations(0),
  _inner_iterations(0),
  _system (s),
  _solve_result(INVALID_SOLVE_RESULT)
{
}



std::unique_ptr<DiffSolver> DiffSolver::build (sys_type & s)
{
  return libmesh_make_unique<NewtonSolver>(s);
}



void DiffSolver::reinit ()
{
  // Reset the max_step_size and max_residual_norm for a new mesh
  max_solution_norm = 0.;
  max_residual_norm = 0.;
}



void DiffSolver::init ()
{
  // Reset the max_step_size and max_residual_norm for a new problem
  max_solution_norm = 0.;
  max_residual_norm = 0.;
}

} // namespace libMesh
