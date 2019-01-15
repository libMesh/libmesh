// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/linear_solver.h"
#include "libmesh/time_solver.h"
#include "libmesh/no_solution_history.h"

namespace libMesh
{



TimeSolver::TimeSolver (sys_type & s)
  : quiet (true),
    reduce_deltat_on_diffsolver_failure (0),
    _diff_solver (),
    _linear_solver (),
    _system (s),
    solution_history(new NoSolutionHistory()), // Default setting for solution_history
    _is_adjoint (false)
{
}



TimeSolver::~TimeSolver ()
{
}



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

void TimeSolver::advance_timestep ()
{
}

void TimeSolver::adjoint_advance_timestep ()
{
}

void TimeSolver::retrieve_timestep ()
{
}

} // namespace libMesh
