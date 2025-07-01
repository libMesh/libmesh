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
 * This class implements the Smoothness indicator.
 *
 * \author ---
 * \author ---
 * \date 2025
 */
class SmoothnessEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Defaults to H1 seminorm.  All Hilbert norms and
   * seminorms should be supported now.  W1,p and W2,p norms would
   * be natural to support if any contributors make the effort.
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
   * This function uses the Patch Recovery error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
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
