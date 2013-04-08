// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_order.h"
#include "libmesh/patch.h"
#include "libmesh/point.h"
#include "libmesh/elem_range.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

// Forward Declarations
class Elem;


/**
 * This class implements the Patch Recovery error indicator.
 *
 *
 * @author Varis Carey, Benjamin S. Kirk, 2004.
 */
class PatchRecoveryErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Defaults to H1 seminorm.  All Hilbert norms and
   * seminorms should be supported now.  W1,p and W2,p norms would
   * be natural to support if any contributors make the effort.
   */
  PatchRecoveryErrorEstimator(const Parallel::Communicator &comm = libMesh::CommWorld) :
    ErrorEstimator(comm),
    target_patch_size(20),
    patch_growth_strategy(&Patch::add_local_face_neighbors),
    patch_reuse(true)
  { error_norm = H1_SEMINORM; }

  /**
   * Destructor.
   */
  ~PatchRecoveryErrorEstimator() {}


  /**
   * This function uses the Patch Recovery error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       const NumericVector<Number>* solution_vector = NULL,
			       bool estimate_parent_error = false);

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

  void set_patch_reuse (bool );

 protected:

  /**
   * Returns the spectral polynomial basis function values at a point x,y,z
   */

  static std::vector<Real> specpoly(const unsigned int dim,
				    const Order order,
				    const Point p,
				    const unsigned int matsize);

  bool patch_reuse ;

 private:

  /**
   * Class to compute the error contribution for a range
   * of elements. May be executed in parallel on separate threads.
   */
  class EstimateError
  {
  public:
    EstimateError (const System& sys,
		   const PatchRecoveryErrorEstimator &ee,
		   ErrorVector& epc) :
      system(sys),
      error_estimator(ee),
      error_per_cell(epc)
    {}

    void operator()(const ConstElemRange &range) const;

    /**
     * Function to set the boolean patch_reuse in case the user
     * wants to change the default behaviour of patch_recovery_error_estimator
     */

  private:

    const System &system;
    const PatchRecoveryErrorEstimator &error_estimator;
    ErrorVector &error_per_cell;
  };

  friend class EstimateError;
};


} // namespace libMesh


#endif // LIBMESH_PATCH_RECOVERY_ERROR_ESTIMATOR_H
