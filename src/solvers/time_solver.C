// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/time_solver.h"

#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/linear_solver.h"
#include "libmesh/no_solution_history.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/error_vector.h"

namespace libMesh
{

TimeSolver::TimeSolver (sys_type & s)
  : quiet (true),
    reduce_deltat_on_diffsolver_failure (0),
    _diff_solver (),
    _linear_solver (),
    _system (s),
    solution_history(libmesh_make_unique<NoSolutionHistory>()),
    last_deltat (s.deltat),
    _is_adjoint (false)
{
}



TimeSolver::~TimeSolver () = default;



void TimeSolver::reinit ()
{
  libmesh_assert(this->diff_solver().get());
  libmesh_assert_equal_to (&(this->diff_solver()->system()), &(this->system()));
  this->diff_solver()->reinit();

  libmesh_assert(this->linear_solver().get());
  this->linear_solver()->clear();
  if (libMesh::on_command_line("--solver-system-names"))
    this->linear_solver()->init((_system.name()+"_").c_str());
  else
    this->linear_solver()->init();

  this->_linear_solver->init_names(_system);
}



void TimeSolver::init ()
{
  // If the user hasn't given us a solver to use,
  // just build a default solver
  if (this->diff_solver().get() == nullptr)
    this->diff_solver() = DiffSolver::build(_system);

  if (this->linear_solver().get() == nullptr)
    this->linear_solver() = LinearSolver<Number>::build(_system.comm());
}



void TimeSolver::init_data ()
{
  this->diff_solver()->init();

  if (libMesh::on_command_line("--solver-system-names"))
    this->linear_solver()->init((_system.name()+"_").c_str());
  else
    this->linear_solver()->init();

  this->linear_solver()->init_names(_system);
}



void TimeSolver::solve ()
{
  libmesh_assert(this->diff_solver().get());
  libmesh_assert_equal_to (&(this->diff_solver()->system()), &(this->system()));
  this->diff_solver()->solve();
}


void TimeSolver::set_solution_history (const SolutionHistory & _solution_history)
{
  solution_history = _solution_history.clone();
}

SolutionHistory & TimeSolver::get_solution_history ()
{
  return *solution_history;
}

void TimeSolver::advance_timestep ()
{
}

std::pair<unsigned int, Real> TimeSolver::adjoint_solve (const QoISet & qoi_indices)
{
  libmesh_assert(this->diff_solver().get());
  libmesh_assert_equal_to (&(this->diff_solver()->system()), &(this->system()));

  return this->_system.ImplicitSystem::adjoint_solve(qoi_indices);
}

void TimeSolver::integrate_qoi_timestep()
{
  libmesh_not_implemented();
}

void TimeSolver::integrate_adjoint_sensitivity(const QoISet & /* qois */, const ParameterVector & /* parameter_vector */, SensitivityData & /* sensitivities */)
{
  libmesh_not_implemented();
}

#ifdef LIBMESH_ENABLE_AMR
void TimeSolver::integrate_adjoint_refinement_error_estimate
  (AdjointRefinementEstimator & /* adjoint_refinement_error_estimator */,
   ErrorVector & /* QoI_elementwise_error */)
{
  libmesh_not_implemented();
}
#endif // LIBMESH_ENABLE_AMR

Real TimeSolver::last_completed_timestep_size()
{
  return last_deltat;
}

void TimeSolver::adjoint_advance_timestep ()
{
}

void TimeSolver::retrieve_timestep ()
{
}

} // namespace libMesh
