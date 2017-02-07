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


// C++ includes
#include <iostream>
#include <iomanip>
#include <sstream>

// Local Includes
#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/error_vector.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/system_norm.h"
#include "libmesh/qoi_set.h"


namespace libMesh
{

//-----------------------------------------------------------------
// AdjointResidualErrorEstimator implementations
AdjointResidualErrorEstimator::AdjointResidualErrorEstimator () :
  ErrorEstimator(),
  error_plot_suffix(),
  _primal_error_estimator(new PatchRecoveryErrorEstimator()),
  _dual_error_estimator(new PatchRecoveryErrorEstimator()),
  _qoi_set(QoISet())
{
}



void AdjointResidualErrorEstimator::estimate_error (const System & _system,
                                                    ErrorVector & error_per_cell,
                                                    const NumericVector<Number> * solution_vector,
                                                    bool estimate_parent_error)
{
  LOG_SCOPE("estimate_error()", "AdjointResidualErrorEstimator");

  // The current mesh
  const MeshBase & mesh = _system.get_mesh();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.max_elem_id());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Get the number of variables in the system
  unsigned int n_vars = _system.n_vars();

  // We need to make a map of the pointer to the solution vector
  std::map<const System *, const NumericVector<Number> *>solutionvecs;
  solutionvecs[&_system] = _system.solution.get();

  // Solve the dual problem if we have to
  if (!_system.is_adjoint_already_solved())
    {
      // FIXME - we'll need to change a lot of APIs to make this trick
      // work with a const System...
      System &  system = const_cast<System &>(_system);
      system.adjoint_solve(_qoi_set);
    }

  // Flag to check whether we have not been asked to weight the variable error contributions in any specific manner
  bool error_norm_is_identity = error_norm.is_identity();

  // Create an ErrorMap/ErrorVector to store the primal, dual and total_dual variable errors
  ErrorMap primal_errors_per_cell;
  ErrorMap dual_errors_per_cell;
  ErrorMap total_dual_errors_per_cell;
  // Allocate ErrorVectors to this map if we're going to use it
  if (!error_norm_is_identity)
    for(unsigned int v = 0; v < n_vars; v++)
      {
        primal_errors_per_cell[std::make_pair(&_system, v)] = new ErrorVector;
        dual_errors_per_cell[std::make_pair(&_system, v)] = new ErrorVector;
        total_dual_errors_per_cell[std::make_pair(&_system, v)] = new ErrorVector;
      }
  ErrorVector primal_error_per_cell;
  ErrorVector dual_error_per_cell;
  ErrorVector total_dual_error_per_cell;

  // Have we been asked to weight the variable error contributions in any specific manner
  if(!error_norm_is_identity) // If we do
    {
      // Estimate the primal problem error for each variable
      _primal_error_estimator->estimate_errors
        (_system.get_equation_systems(), primal_errors_per_cell, &solutionvecs, estimate_parent_error);
    }
  else // If not
    {
      // Just get the combined error estimate
      _primal_error_estimator->estimate_error
        (_system, primal_error_per_cell, solution_vector, estimate_parent_error);
    }

  // Sum and weight the dual error estimate based on our QoISet
  for (std::size_t i = 0; i != _system.qoi.size(); ++i)
    {
      if (_qoi_set.has_index(i))
        {
          // Get the weight for the current QoI
          Real error_weight = _qoi_set.weight(i);

          // We need to make a map of the pointer to the adjoint solution vector
          std::map<const System *, const NumericVector<Number> *>adjointsolutionvecs;
          adjointsolutionvecs[&_system] = &_system.get_adjoint_solution(i);

          // Have we been asked to weight the variable error contributions in any specific manner
          if(!error_norm_is_identity) // If we have
            {
              _dual_error_estimator->estimate_errors
                (_system.get_equation_systems(), dual_errors_per_cell, &adjointsolutionvecs,
                 estimate_parent_error);
            }
          else // If not
            {
              // Just get the combined error estimate
              _dual_error_estimator->estimate_error
                (_system, dual_error_per_cell, &(_system.get_adjoint_solution(i)), estimate_parent_error);
            }

          std::size_t error_size;

          // Get the size of the first ErrorMap vector; this will give us the number of elements
          if(!error_norm_is_identity) // If in non default weights case
            {
              error_size = dual_errors_per_cell[std::make_pair(&_system, 0)]->size();
            }
          else // If in the standard default weights case
            {
              error_size = dual_error_per_cell.size();
            }

          // Resize the ErrorVector(s)
          if(!error_norm_is_identity)
            {
              // Loop over variables
              for(unsigned int v = 0; v < n_vars; v++)
                {
                  libmesh_assert(!total_dual_errors_per_cell[std::make_pair(&_system, v)]->size() ||
                                 total_dual_errors_per_cell[std::make_pair(&_system, v)]->size() == error_size) ;
                  total_dual_errors_per_cell[std::make_pair(&_system, v)]->resize(error_size);
                }
            }
          else
            {
              libmesh_assert(!total_dual_error_per_cell.size() ||
                             total_dual_error_per_cell.size() == error_size);
              total_dual_error_per_cell.resize(error_size);
            }

          for (std::size_t e = 0; e != error_size; ++e)
            {
              // Have we been asked to weight the variable error contributions in any specific manner
              if(!error_norm_is_identity) // If we have
                {
                  // Loop over variables
                  for(unsigned int v = 0; v < n_vars; v++)
                    {
                      // Now fill in total_dual_error ErrorMap with the weight
                      (*total_dual_errors_per_cell[std::make_pair(&_system, v)])[e] +=
                        static_cast<ErrorVectorReal>
                        (error_weight *
                         (*dual_errors_per_cell[std::make_pair(&_system, v)])[e]);
                    }
                }
              else // If not
                {
                  total_dual_error_per_cell[e] +=
                    static_cast<ErrorVectorReal>(error_weight * dual_error_per_cell[e]);
                }
            }
        }
    }

  // Do some debugging plots if requested
  if (!error_plot_suffix.empty())
    {
      if(!error_norm_is_identity) // If we have
        {
          // Loop over variables
          for(unsigned int v = 0; v < n_vars; v++)
            {
              std::ostringstream primal_out;
              std::ostringstream dual_out;
              primal_out << "primal_" << error_plot_suffix << ".";
              dual_out << "dual_" << error_plot_suffix << ".";

              primal_out << std::setw(1)
                         << std::setprecision(0)
                         << std::setfill('0')
                         << std::right
                         << v;

              dual_out << std::setw(1)
                       << std::setprecision(0)
                       << std::setfill('0')
                       << std::right
                       << v;

              (*primal_errors_per_cell[std::make_pair(&_system, v)]).plot_error(primal_out.str(), _system.get_mesh());
              (*total_dual_errors_per_cell[std::make_pair(&_system, v)]).plot_error(dual_out.str(), _system.get_mesh());

              primal_out.clear();
              dual_out.clear();
            }
        }
      else // If not
        {
          std::ostringstream primal_out;
          std::ostringstream dual_out;
          primal_out << "primal_" << error_plot_suffix ;
          dual_out << "dual_" << error_plot_suffix ;

          primal_error_per_cell.plot_error(primal_out.str(), _system.get_mesh());
          total_dual_error_per_cell.plot_error(dual_out.str(), _system.get_mesh());

          primal_out.clear();
          dual_out.clear();
        }
    }

  // Weight the primal error by the dual error using the system norm object
  // FIXME: we ought to thread this
  for (std::size_t i=0; i != error_per_cell.size(); ++i)
    {
      // Have we been asked to weight the variable error contributions in any specific manner
      if(!error_norm_is_identity) // If we do
        {
          // Create Error Vectors to pass to calculate_norm
          std::vector<Real> cell_primal_error;
          std::vector<Real> cell_dual_error;

          for(unsigned int v = 0; v < n_vars; v++)
            {
              cell_primal_error.push_back((*primal_errors_per_cell[std::make_pair(&_system, v)])[i]);
              cell_dual_error.push_back((*total_dual_errors_per_cell[std::make_pair(&_system, v)])[i]);
            }

          error_per_cell[i] =
            static_cast<ErrorVectorReal>
            (error_norm.calculate_norm(cell_primal_error, cell_dual_error));
        }
      else // If not
        {
          error_per_cell[i] = primal_error_per_cell[i]*total_dual_error_per_cell[i];
        }
    }

  // Deallocate the ErrorMap contents if we allocated them earlier
  if (!error_norm_is_identity)
    for(unsigned int v = 0; v < n_vars; v++)
      {
        delete primal_errors_per_cell[std::make_pair(&_system, v)];
        delete dual_errors_per_cell[std::make_pair(&_system, v)];
        delete total_dual_errors_per_cell[std::make_pair(&_system, v)];
      }
}

} // namespace libMesh
