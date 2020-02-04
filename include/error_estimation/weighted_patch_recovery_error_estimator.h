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



#ifndef LIBMESH_WEIGHTED_PATCH_RECOVERY_ERROR_ESTIMATOR_H
#define LIBMESH_WEIGHTED_PATCH_RECOVERY_ERROR_ESTIMATOR_H

// C++ includes
#include <vector>

// Local Includes
#include "libmesh/elem_range.h"
#include "libmesh/error_estimator.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/patch.h"


namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;


/**
 * This class implements the Patch Recovery error indicator.
 *
 * \author Vikram Garg
 * \date 2012
 */
class WeightedPatchRecoveryErrorEstimator : public PatchRecoveryErrorEstimator
{
public:

  /**
   * Constructor.  Defaults to H1 seminorm.  All Hilbert norms and
   * seminorms should be supported now.  W1,p and W2,p norms would
   * be natural to support if any contributors make the effort.
   */
  WeightedPatchRecoveryErrorEstimator() = default;

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  WeightedPatchRecoveryErrorEstimator (const WeightedPatchRecoveryErrorEstimator &) = default;
  WeightedPatchRecoveryErrorEstimator (WeightedPatchRecoveryErrorEstimator &&) = default;
  WeightedPatchRecoveryErrorEstimator & operator= (const WeightedPatchRecoveryErrorEstimator &) = default;
  WeightedPatchRecoveryErrorEstimator & operator= (WeightedPatchRecoveryErrorEstimator &&) = default;
  virtual ~WeightedPatchRecoveryErrorEstimator() = default;

  /**
   * This function uses the Patch Recovery error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = nullptr,
                               bool estimate_parent_error = false) override;

  /**
   * Vector of fem function base pointers, the user will fill this in
   * with pointers to the appropriate weight functions.
   */
  std::vector<FEMFunctionBase<Number> *> weight_functions;

  virtual ErrorEstimatorType type() const override;

private:

  /**
   * Class to compute the error contribution for a range
   * of elements. May be executed in parallel on separate threads.
   */
  class EstimateError
  {
  public:
    EstimateError (const System & sys,
                   const WeightedPatchRecoveryErrorEstimator & ee,
                   ErrorVector & epc) :
      system(sys),
      error_estimator(ee),
      error_per_cell(epc)
    {}

    void operator()(const ConstElemRange & range) const;

    /**
     * Function to set the boolean patch_reuse in case the user
     * wants to change the default behaviour of patch_recovery_error_estimator
     */

  private:

    const System & system;
    const WeightedPatchRecoveryErrorEstimator & error_estimator;
    ErrorVector & error_per_cell;
  };

  friend class EstimateError;
};


} // namespace libMesh


#endif // LIBMESH_WEIGHTED_PATCH_RECOVERY_ERROR_ESTIMATOR_H
