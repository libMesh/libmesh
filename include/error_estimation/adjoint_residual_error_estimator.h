// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ADJOINT_RESIDUAL_ERROR_ESTIMATOR_H
#define LIBMESH_ADJOINT_RESIDUAL_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/error_estimator.h"
#include "libmesh/qoi_set.h"

// C++ includes
#include <cstddef>
#include <string>
#include <vector>

// Forward Declarations




namespace libMesh
{


/**
 * This class implements a goal oriented error indicator, by weighting
 * residual-based estimates from the primal problem against estimates
 * from the adjoint problem.
 *
 * This is based on a trick suggested by Brian Carnes, (first proposed by
 * Babuska and Miller in 1984) but if it
 * doesn't actually work then the misunderstanding or
 * misimplementation will be the fault of Roy Stogner.  It's also
 * Roy's fault there's no literature reference here yet.
 *
 * \author Roy H. Stogner
 * \date 2009
 */
class AdjointResidualErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for picking default subestimators.
   */
  AdjointResidualErrorEstimator();

  /**
   * Destructor.
   */
  ~AdjointResidualErrorEstimator() {}

  /**
   * Access to the "subestimator" (default: PatchRecovery) to use on
   * the primal/forward solution
   */
  UniquePtr<ErrorEstimator> & primal_error_estimator() { return _primal_error_estimator; }

  /**
   * Access to the "subestimator" (default: PatchRecovery) to use on
   * the dual/adjoint solution
   */
  UniquePtr<ErrorEstimator> & dual_error_estimator() { return _dual_error_estimator; }

  /**
   * Access to the QoISet (default: weight all QoIs equally) to use
   * when computing errors
   */
  QoISet & qoi_set() { return _qoi_set; }

  /**
   * Access to the QoISet (default: weight all QoIs equally) to use
   * when computing errors
   */
  const QoISet & qoi_set() const { return _qoi_set; }

  /**
   * To aid in investigating error estimator behavior, set this string
   * to a suffix with which to plot (prefixed by "primal_" or "dual_")
   * the subestimator results.  The suffix should end with a file
   * extension (e.g. ".gmv") that the ErrorVector::plot_error
   * recognizes.
   */
  std::string error_plot_suffix;

  /**
   * Compute the adjoint-weighted error on each element and place it
   * in the \p error_per_cell vector.  Note that this->error_norm is
   * ignored; the error estimate is in the seminorm given by the
   * absolute value of the error in the quantity of interest
   * functional.  The primal and dual subestimator error_norm values
   * are used, and should be chosen appropriately for your model.
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = libmesh_nullptr,
                               bool estimate_parent_error = false) libmesh_override;

  virtual ErrorEstimatorType type() const libmesh_override
  { return ADJOINT_RESIDUAL;}

protected:

  /**
   * An error estimator for the forward problem
   */
  UniquePtr<ErrorEstimator> _primal_error_estimator;

  /**
   * An error estimator for the adjoint problem
   */
  UniquePtr<ErrorEstimator> _dual_error_estimator;

  /**
   * A QoISet to handle cases with multiple QoIs available
   */
  QoISet _qoi_set;
};


} // namespace libMesh

#endif // LIBMESH_ADJOINT_RESIDUAL_ERROR_ESTIMATOR_H
