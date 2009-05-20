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


// C++ includes


// Local Includes
#include "adjoint_residual_error_estimator.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "system.h"


//-----------------------------------------------------------------
// AdjointResidualErrorEstimator implementations
AdjointResidualErrorEstimator::AdjointResidualErrorEstimator () :
  adjoint_already_solved(false),
  _primal_error_estimator(new KellyErrorEstimator()),
  _dual_error_estimator(new KellyErrorEstimator())
{
}



void AdjointResidualErrorEstimator::estimate_error (const System& _system,
						    ErrorVector& error_per_cell,
					            bool estimate_parent_error)
{
  START_LOG("estimate_error()", "AdjointResidualErrorEstimator");

  // FIXME - we'll need to change a lot of APIs to make this trick
  // work with a const System...
  System&  system = const_cast<System&>(_system);

  // Start by estimating the primal problem error
  _primal_error_estimator->estimate_error
    (system, error_per_cell, estimate_parent_error);

  // Solve the dual problem if we have to
  if (!adjoint_already_solved)
    system.adjoint_solve();

  // Swap the system solution and dual solution, so our next error
  // estimate will use the latter
  system.solution->swap(system.get_adjoint_solution());

  // Get a separate estimate of the dual problem error
  ErrorVector dual_error_per_cell;
  _dual_error_estimator->estimate_error
    (system, dual_error_per_cell, estimate_parent_error);

  // Swap the system solution and dual solution back to their proper
  // places
  system.solution->swap(system.get_adjoint_solution());

  // Weight the primal error by the dual error
  // FIXME: we ought to thread this
  for (unsigned int i=0; i != error_per_cell.size(); ++i)
    error_per_cell[i] *= dual_error_per_cell[i];

  STOP_LOG("estimate_error()", "AdjointResidualErrorEstimator");
}
