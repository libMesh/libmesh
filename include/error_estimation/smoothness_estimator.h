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



#ifndef LIBMESH_SMOOTHNESS_ESTIMATOR_H
#define LIBMESH_SMOOTHNESS_ESTIMATOR_H

// Local Includes
#include "libmesh/error_estimator.h"
#include "libmesh/point.h"
#include "libmesh/elem_range.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

// Forward Declarations
class Elem;
enum Order : int;

/**
 * This class implements the Smoothness estimate.
 *
 * \author Arjun Kulathuvayal
 * \author Varis Carey
 * \author Benjamin S. Kirk
 * \date 2025
 */
class SmoothnessEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.
   */
  SmoothnessEstimator();

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  SmoothnessEstimator (const SmoothnessEstimator &) = default;
  SmoothnessEstimator (SmoothnessEstimator &&) = default;
  SmoothnessEstimator & operator= (const SmoothnessEstimator &) = default;
  SmoothnessEstimator & operator= (SmoothnessEstimator &&) = default;
  virtual ~SmoothnessEstimator() = default;

  /**
   * This function uses the Legendre expansion of solution to estimate coefficient
   * decay to quantify the solution smoothness on each cell.
   * The estimated smoothness is output in the vector (stored as ErrorVector)
   * For element order 1, the least square fit of log|order| vs log |coefficients|
   * fails. This leads to pure h refinement.
   * \p smoothness_per_cell
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & smoothness_per_cell,
                               const NumericVector<Number> * solution_vector = nullptr,
                               bool estimate_parent_error = false) override;

  /**
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.  The default \p extraorder is 1, because properly
   * integrating L2 error requires integrating the squares of terms
   * with order p+1, and 2p+2 is 1 higher than what we default to
   * using for reasonable mass matrix integration.
   */
  void extra_quadrature_order (const int extraorder)
  { _extra_order = extraorder; }

  virtual ErrorEstimatorType type() const override;

protected:

  /**
   * \returns The Legendre polynomial basis function values at a point (x,y,z).
   */
  static std::vector<Real> legepoly(const unsigned int dim,
                                    const Order order,
                                    const Point p,
                                    const unsigned int matsize);

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;

  /**
   * Computes slop in a linear regression
   */
  static Real compute_slope(int N, Real Sx, Real Sy, Real Sxx, Real Sxy);
  
private:

  /**
   * Class to compute the error contribution for a range
   * of elements. May be executed in parallel on separate threads.
   */
  class EstimateError
  {
  public:
    EstimateError (const System & sys,
                   const SmoothnessEstimator & ee,
                   ErrorVector & epc) :
      system(sys),
      error_estimator(ee),
      smoothness_per_cell(epc)
    {}

    void operator()(const ConstElemRange & range) const;

  private:
    const System & system;
    const SmoothnessEstimator & error_estimator;
    ErrorVector & smoothness_per_cell;
  };

  friend class EstimateError;
};


} // namespace libMesh


#endif // LIBMESH_SMOOTHNESS_ESTIMATOR_H
