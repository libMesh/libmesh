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



#ifndef LIBMESH_PATCH_RECOVERY_ERROR_ESTIMATOR_H
#define LIBMESH_PATCH_RECOVERY_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/error_estimator.h"
#include "libmesh/patch.h"
#include "libmesh/point.h"
#include "libmesh/elem_range.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum Order : int;
}
#else
#include "libmesh/enum_order.h"
#endif

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;

/**
 * This class implements the Patch Recovery error indicator.
 *
 * \author Varis Carey
 * \author Benjamin S. Kirk
 * \date 2004
 */
class PatchRecoveryErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Defaults to H1 seminorm.  All Hilbert norms and
   * seminorms should be supported now.  W1,p and W2,p norms would
   * be natural to support if any contributors make the effort.
   */
  PatchRecoveryErrorEstimator();

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  PatchRecoveryErrorEstimator (const PatchRecoveryErrorEstimator &) = default;
  PatchRecoveryErrorEstimator (PatchRecoveryErrorEstimator &&) = default;
  PatchRecoveryErrorEstimator & operator= (const PatchRecoveryErrorEstimator &) = default;
  PatchRecoveryErrorEstimator & operator= (PatchRecoveryErrorEstimator &&) = default;
  virtual ~PatchRecoveryErrorEstimator() = default;

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
   * The PatchErrorEstimator will build patches of at least this many
   * elements to perform estimates
   */
  unsigned int target_patch_size;

  /**
   * The PatchErrorEstimator will use this pointer to a Patch member
   * function when growing patches.  The default strategy used is
   * Patch::add_local_face_neighbors.
   * Patch::add_local_point_neighbors may be more reliable but slower.
   */
  Patch::PMF patch_growth_strategy;

  void set_patch_reuse (bool);

  virtual ErrorEstimatorType type() const override;

protected:

  /**
   * \returns The spectral polynomial basis function values at a point (x,y,z).
   */
  static std::vector<Real> specpoly(const unsigned int dim,
                                    const Order order,
                                    const Point p,
                                    const unsigned int matsize);

  bool patch_reuse;

private:

  /**
   * Class to compute the error contribution for a range
   * of elements. May be executed in parallel on separate threads.
   */
  class EstimateError
  {
  public:
    EstimateError (const System & sys,
                   const PatchRecoveryErrorEstimator & ee,
                   ErrorVector & epc) :
      system(sys),
      error_estimator(ee),
      error_per_cell(epc)
    {}

    void operator()(const ConstElemRange & range) const;

  private:
    const System & system;
    const PatchRecoveryErrorEstimator & error_estimator;
    ErrorVector & error_per_cell;
  };

  friend class EstimateError;
};


} // namespace libMesh


#endif // LIBMESH_PATCH_RECOVERY_ERROR_ESTIMATOR_H
