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
#include "patch_recovery_error_estimator.h"
#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "system.h"
#include "qoi_set.h"


//-----------------------------------------------------------------
// AdjointResidualErrorEstimator implementations
AdjointResidualErrorEstimator::AdjointResidualErrorEstimator () :
  adjoint_already_solved(false),
  error_plot_suffix(),
  _primal_error_estimator(new PatchRecoveryErrorEstimator()),
  _dual_error_estimator(new PatchRecoveryErrorEstimator()),
  _qoi_set(QoISet())
{
}



void AdjointResidualErrorEstimator::estimate_error (const System& _system,
						    ErrorVector& error_per_cell,
						    const NumericVector<Number>* solution_vector,
					            bool estimate_parent_error)
{
  START_LOG("estimate_error()", "AdjointResidualErrorEstimator");

  // Start by estimating the primal problem error
  _primal_error_estimator->estimate_error
    (_system, error_per_cell, solution_vector, estimate_parent_error);

  // Solve the dual problem if we have to
  if (!adjoint_already_solved)
    {
      // FIXME - we'll need to change a lot of APIs to make this trick
      // work with a const System...
      System&  system = const_cast<System&>(_system);
      system.adjoint_solve(_qoi_set);
    }

  // This bookkeeping should now be taken care of in subestimators
  // when they see a non-default solution_vector
  //
  // // Swap the system solution and dual solution, so our next error
  // // estimate will use the latter
  // // system.solution->swap(system.get_adjoint_solution());

  // // Don't forget to keep the current_local_solution consistent, as
  // // error estimators are likely to use that!
  // // system.update();

  // Get a separate estimate of the dual problem error
  ErrorVector total_dual_error_per_cell;

  // Sum and weight this estimate based on our QoISet
  for (unsigned int i = 0; i != _system.qoi.size(); ++i)
    if (_qoi_set.has_index(i))
      {
        Real error_weight = _qoi_set.weight(i);

        ErrorVector dual_error_per_cell;
        _dual_error_estimator->estimate_error
          (_system, dual_error_per_cell, &(_system.get_adjoint_solution(i)),
           estimate_parent_error);

        unsigned int error_size = dual_error_per_cell.size();
        libmesh_assert(!total_dual_error_per_cell.size() ||
                       total_dual_error_per_cell.size() == error_size);
        total_dual_error_per_cell.resize(error_size);
        for (unsigned int e = 0; e != error_size; ++e)
          total_dual_error_per_cell[e] += 
            error_weight * dual_error_per_cell[e];
      }

  // Do some debugging plots if requested
  if (!error_plot_suffix.empty())
    {
      std::string primal_out = "primal_";
      std::string dual_out = "dual_";
      primal_out += error_plot_suffix;
      dual_out += error_plot_suffix;
      error_per_cell.plot_error(primal_out, _system.get_mesh());
      total_dual_error_per_cell.plot_error(dual_out, _system.get_mesh());
    }

  // More bookkeeping for subestimators to do
  //
  // // Swap the system solution and dual solution back to their proper
  // // places
  // system.solution->swap(system.get_adjoint_solution());

  // // Don't forget to fix the current_local_solution
  // system.update();

  // Weight the primal error by the dual error
  // FIXME: we ought to thread this
  for (unsigned int i=0; i != error_per_cell.size(); ++i)
    error_per_cell[i] *= total_dual_error_per_cell[i];

  STOP_LOG("estimate_error()", "AdjointResidualErrorEstimator");
}
