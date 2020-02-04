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



#ifndef LIBMESH_EXACT_ERROR_ESTIMATOR_H
#define LIBMESH_EXACT_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/error_estimator.h"
#include "libmesh/function_base.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/fe_base_forward.h"

// C++ includes
#include <cstddef>
#include <string>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
typedef FEGenericBase<Real> FEBase;
class MeshFunction;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;
class Parameters;
template <typename T> class DenseVector;

/**
 * This class implements an "error estimator"
 * based on the difference between the approximate
 * and exact solution.  In theory the quadrature error
 * in this estimate should be much lower than the
 * approximation error in other estimates, so this
 * estimator can be used to calculate effectivity.
 *
 * \author Roy Stogner
 * \date 2006
 */
class ExactErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to nullptr, and defaulting the norm type to H1.
   */
  ExactErrorEstimator();

  /**
   * This class cannot be (default) copy constructed/assigned because
   * it has containers of unique_ptrs. Explicitly deleting these functions is
   * the best way to document this fact.
   */
  ExactErrorEstimator (const ExactErrorEstimator &) = delete;
  ExactErrorEstimator & operator= (const ExactErrorEstimator &) = delete;

  /**
   * Defaulted move ctor, move assignment operator, and destructor.
   */
  ExactErrorEstimator (ExactErrorEstimator &&) = default;
  ExactErrorEstimator & operator= (ExactErrorEstimator &&) = default;
  virtual ~ExactErrorEstimator() = default;

  /**
   * Clone and attach arbitrary functors which compute the exact
   * values of the EquationSystems' solutions at any point.
   */
  void attach_exact_values (std::vector<FunctionBase<Number> *> f);

  /**
   * Clone and attach an arbitrary functor which computes the exact
   * value of the system \p sys_num solution at any point.
   */
  void attach_exact_value (unsigned int sys_num,
                           FunctionBase<Number> * f);

  /**
   * Attach an arbitrary function which computes the exact value of
   * the solution at any point.
   */
  typedef Number (*ValueFunctionPointer)(const Point & p,
                                         const Parameters & Parameters,
                                         const std::string & sys_name,
                                         const std::string & unknown_name);
  void attach_exact_value (ValueFunctionPointer fptr);

  /**
   * Clone and attach arbitrary functors which compute the exact
   * gradients of the EquationSystems' solutions at any point.
   */
  void attach_exact_derivs (std::vector<FunctionBase<Gradient> *> g);

  /**
   * Clone and attach an arbitrary functor which computes the exact
   * gradient of the system \p sys_num solution at any point.
   */
  void attach_exact_deriv (unsigned int sys_num,
                           FunctionBase<Gradient> * g);

  /**
   * Attach an arbitrary function which computes the exact gradient of
   * the solution at any point.
   */
  typedef Gradient (*GradientFunctionPointer)(const Point & p,
                                              const Parameters & parameters,
                                              const std::string & sys_name,
                                              const std::string & unknown_name);
  void attach_exact_deriv (GradientFunctionPointer gptr);

  /**
   * Clone and attach arbitrary functors which compute the exact
   * second derivatives of the EquationSystems' solutions at any point.
   */
  void attach_exact_hessians (std::vector<FunctionBase<Tensor> *> h);

  /**
   * Clone and attach an arbitrary functor which computes the exact
   * second derivatives of the system \p sys_num solution at any point.
   */
  void attach_exact_hessian (unsigned int sys_num,
                             FunctionBase<Tensor> * h);

  /**
   * Attach an arbitrary function which computes the exact second
   * derivatives of the solution at any point.
   */
  typedef Tensor (*HessianFunctionPointer)(const Point & p,
                                           const Parameters & parameters,
                                           const std::string & sys_name,
                                           const std::string & unknown_name);
  void attach_exact_hessian (HessianFunctionPointer hptr);

  /**
   * Attach function similar to system.h which
   * allows the user to attach a second EquationSystems
   * object with a reference fine grid solution.
   */
  void attach_reference_solution (EquationSystems * es_fine);


  /**
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.
   */
  void extra_quadrature_order (const int extraorder)
  { _extra_order = extraorder; }


  // Bring the base class functionality into the name lookup
  // procedure.  This allows for alternative calling formats
  // defined in the base class.  Thanks Wolfgang.
  // GCC 2.95.3 cannot compile such code.  Since it was not really
  // essential to the functioning of this class, it's been removed.
  // using ErrorEstimator::estimate_error;

  /**
   * This function uses the exact solution function
   * to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = nullptr,
                               bool estimate_parent_error = false) override;

  virtual ErrorEstimatorType type() const override;

private:

  /**
   * Function pointer to user-provided function which
   * computes the exact value of the solution.
   */
  ValueFunctionPointer _exact_value;

  /**
   * Function pointer to user-provided function which
   * computes the exact derivative of the solution.
   */
  GradientFunctionPointer _exact_deriv;

  /**
   * Function pointer to user-provided function which
   * computes the exact hessian of the solution.
   */
  HessianFunctionPointer _exact_hessian;

  /**
   * User-provided functors which compute the exact value of the
   * solution for each system.
   */
  std::vector<std::unique_ptr<FunctionBase<Number>>> _exact_values;

  /**
   * User-provided functors which compute the exact derivative of the
   * solution for each system.
   */
  std::vector<std::unique_ptr<FunctionBase<Gradient>>> _exact_derivs;

  /**
   * User-provided functors which compute the exact hessians of the
   * solution for each system.
   */
  std::vector<std::unique_ptr<FunctionBase<Tensor>>> _exact_hessians;

  /**
   * Constant pointer to the \p EquationSystems object
   * containing a fine grid solution.
   */
  EquationSystems * _equation_systems_fine;

  /**
   * Helper method for calculating on each element
   */
  Real find_squared_element_error (const System & system,
                                   const std::string & var_name,
                                   const Elem * elem,
                                   const DenseVector<Number> & Uelem,
                                   FEBase * fe,
                                   MeshFunction * fine_values) const;

  /**
   * Helper method for cleanup
   */
  void clear_functors ();

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;
};


} // namespace libMesh

#endif // LIBMESH_EXACT_ERROR_ESTIMATOR_H
