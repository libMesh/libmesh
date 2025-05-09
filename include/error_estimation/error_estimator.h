// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ERROR_ESTIMATOR_H
#define LIBMESH_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/system_norm.h"

// C++ includes
#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace libMesh
{

// Forward Declarations
class ErrorVector;
class EquationSystems;
class System;
template <typename T> class NumericVector;

namespace Parallel {
  class Communicator;
}

enum ErrorEstimatorType : int;

/**
 * This class holds functions that will estimate the error
 * in a finite element solution on a given mesh.  These error
 * estimates can be useful in their own right, or may be used
 * to guide adaptive mesh refinement.
 *
 * \note The computed errors are stored as floats rather
 * than doubles since the required precision is low.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
class ErrorEstimator
{
public:

  /**
   * Constructor.  Empty.  Derived classes should reset error_norm as
   * appropriate.
   */
  ErrorEstimator() = default;

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  ErrorEstimator (const ErrorEstimator &) = default;
  ErrorEstimator (ErrorEstimator &&) = default;
  ErrorEstimator & operator= (const ErrorEstimator &) = default;
  ErrorEstimator & operator= (ErrorEstimator &&) = default;
  virtual ~ErrorEstimator() = default;


  /**
   * This pure virtual function must be redefined
   * in derived classes to compute the error for each
   * active element and place it in the "error_per_cell" vector.
   *
   * If solution_vector is not nullptr, the estimator will
   * (if able) attempt to estimate an error in that field
   * instead of in system.solution.
   *
   * If estimate_parent_error is not false, the estimator will (if
   * able) attempt to give a consistent estimate of errors in parent
   * elements that would be generated by coarsening.
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = nullptr,
                               bool estimate_parent_error = false) = 0;

  /**
   * This virtual function can be redefined
   * in derived classes, but by default computes the sum of
   * the error_per_cell for each system in the equation_systems.
   *
   * Currently this function ignores the error_norm member variable,
   * and uses the function argument error_norms instead.
   *
   * This function is named estimate_errors instead of estimate_error
   * because otherwise C++ can get confused.
   */
  virtual void estimate_errors (const EquationSystems & equation_systems,
                                ErrorVector & error_per_cell,
                                const std::map<const System *, SystemNorm> & error_norms,
                                const std::map<const System *, const NumericVector<Number> *> * solution_vectors = nullptr,
                                bool estimate_parent_error = false);

  /**
   * When calculating many error vectors at once, we need a data structure to
   * hold them all
   */
  typedef std::map<std::pair<const System *, unsigned int>, std::unique_ptr<ErrorVector>> ErrorMap;

  /**
   * This virtual function can be redefined
   * in derived classes, but by default it calls estimate_error
   * repeatedly to calculate the requested error vectors.
   *
   * Currently this function ignores the error_norm.weight() values
   * because it calculates each variable's error individually, unscaled.
   *
   * The user selects which errors get computed by filling a map with error
   * vectors: If errors_per_cell[&system][v] exists, it will be filled with the
   * error values in variable \p v of \p system
   */
  virtual void estimate_errors (const EquationSystems & equation_systems,
                                ErrorMap & errors_per_cell,
                                const std::map<const System *, const NumericVector<Number> *> * solution_vectors = nullptr,
                                bool estimate_parent_error = false);

  /**
   * \returns The type for the ErrorEstimator subclass.
   */
  virtual ErrorEstimatorType type() const = 0;

  /**
   * When estimating the error in a single system, the \p error_norm
   * is used to control the scaling and norm choice for each variable.
   * Not all estimators will support all norm choices.  The default
   * scaling is for all variables to be weighted equally.  The default
   * norm choice depends on the error estimator.
   *
   * Part of this functionality was supported via component_scale and
   * sobolev_order in older libMesh versions, and a small part was
   * supported via component_mask in even older versions.  Hopefully
   * the encapsulation here will allow us to avoid changing this API
   * again.
   */
  SystemNorm error_norm;

protected:

  /**
   * This method takes the local error contributions in
   * \p error_per_cell from each processor and combines
   * them to get the global error vector.
   */
  void reduce_error (std::vector<ErrorVectorReal> & error_per_cell,
                     const Parallel::Communicator & comm) const;
};


} // namespace libMesh

#endif // LIBMESH_ERROR_ESTIMATOR_H
