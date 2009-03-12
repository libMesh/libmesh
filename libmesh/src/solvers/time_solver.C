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
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s)
  : quiet (true),
    reduce_deltat_on_diffsolver_failure (0),
    _diff_solver (NULL),
    _adjoint_solver (NULL),
    _system (s)
{
}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::reinit ()
{
  _diff_solver->reinit();
  _adjoint_solver->reinit();
}



void TimeSolver::init ()
{
  // If the user hasn't given us a solver to use,
  // just build a default solver
  if (_diff_solver.get() == NULL)
    _diff_solver = DiffSolver::build(_system);

  if (_adjoint_solver.get() == NULL)
    _adjoint_solver = DiffSolver::build(_system);

  _diff_solver->init();
  _adjoint_solver->init();
}



void TimeSolver::solve ()
{
  _diff_solver->solve();
}



void TimeSolver::adjoint_solve ()
{
  _adjoint_solver->adjoint_solve();
}



void TimeSolver::advance_timestep ()
{
}



void TimeSolver::adjoint_recede_timestep ()
{
}
