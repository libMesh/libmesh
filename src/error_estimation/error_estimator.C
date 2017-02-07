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


// Local include files
#include "libmesh/libmesh_common.h"
#include "libmesh/error_estimator.h"
#include "libmesh/error_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel.h"




namespace libMesh
{

// ErrorEstimator functions
void ErrorEstimator::reduce_error (std::vector<ErrorVectorReal> & error_per_cell,
                                   const Parallel::Communicator & comm) const
{
  // This function must be run on all processors at once
  // parallel_object_only();

  // Each processor has now computed the error contribuions
  // for its local elements.  We may need to sum the vector to
  // recover the error for each element.

  comm.sum(error_per_cell);
}



void ErrorEstimator::estimate_errors(const EquationSystems & equation_systems,
                                     ErrorVector & error_per_cell,
                                     const std::map<const System *, SystemNorm> & error_norms,
                                     const std::map<const System *, const NumericVector<Number> *> * solution_vectors,
                                     bool estimate_parent_error)
{
  SystemNorm old_error_norm = this->error_norm;

  // Sum the error values from each system
  for (unsigned int s = 0; s != equation_systems.n_systems(); ++s)
    {
      ErrorVector system_error_per_cell;
      const System & sys = equation_systems.get_system(s);
      if (error_norms.find(&sys) == error_norms.end())
        this->error_norm = old_error_norm;
      else
        this->error_norm = error_norms.find(&sys)->second;

      const NumericVector<Number> * solution_vector = libmesh_nullptr;
      if (solution_vectors &&
          solution_vectors->find(&sys) != solution_vectors->end())
        solution_vector = solution_vectors->find(&sys)->second;

      this->estimate_error(sys, system_error_per_cell,
                           solution_vector, estimate_parent_error);

      if (s)
        {
          libmesh_assert_equal_to (error_per_cell.size(), system_error_per_cell.size());
          for (std::size_t i=0; i != error_per_cell.size(); ++i)
            error_per_cell[i] += system_error_per_cell[i];
        }
      else
        error_per_cell = system_error_per_cell;
    }

  // Restore our old state before returning
  this->error_norm = old_error_norm;
}



/**
 * FIXME: This is a default implementation - derived classes should
 * reimplement it for efficiency.
 */
void ErrorEstimator::estimate_errors(const EquationSystems & equation_systems,
                                     ErrorMap & errors_per_cell,
                                     const std::map<const System *, const NumericVector<Number> *> * solution_vectors,
                                     bool estimate_parent_error)
{
  SystemNorm old_error_norm = this->error_norm;

  // Find the requested error values from each system
  for (unsigned int s = 0; s != equation_systems.n_systems(); ++s)
    {
      const System & sys = equation_systems.get_system(s);

      unsigned int n_vars = sys.n_vars();

      for (unsigned int v = 0; v != n_vars; ++v)
        {
          // Only fill in ErrorVectors the user asks for
          if (errors_per_cell.find(std::make_pair(&sys, v)) ==
              errors_per_cell.end())
            continue;

          // Calculate error in only one variable
          std::vector<Real> weights(n_vars, 0.0);
          weights[v] = 1.0;
          this->error_norm =
            SystemNorm(std::vector<FEMNormType>(n_vars, old_error_norm.type(v)),
                       weights);

          const NumericVector<Number> * solution_vector = libmesh_nullptr;
          if (solution_vectors &&
              solution_vectors->find(&sys) != solution_vectors->end())
            solution_vector = solution_vectors->find(&sys)->second;

          this->estimate_error
            (sys, *errors_per_cell[std::make_pair(&sys, v)],
             solution_vector, estimate_parent_error);
        }
    }

  // Restore our old state before returning
  this->error_norm = old_error_norm;
}

} // namespace libMesh
