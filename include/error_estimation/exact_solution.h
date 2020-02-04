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

#ifndef LIBMESH_EXACT_SOLUTION_H
#define LIBMESH_EXACT_SOLUTION_H


// Local Includes
#include "libmesh/libmesh_common.h" // for Number

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum FEMNormType : int;
}
#else
#include "libmesh/enum_norm_type.h"
#endif

// C++ includes
#include <map>
#include <vector>
#include <memory>
#include <set>

namespace libMesh
{


// Forward Declarations
template <typename> class PointTempl;
typedef PointTempl<Real> Point;
class EquationSystems;
class Parameters;
class Mesh;
template <typename Output> class FunctionBase;

// Is there any way to simplify this?
// All we need are Tensor and Gradient. - RHS
template <typename T> class TensorValue;
template <typename T> class VectorValue;
typedef TensorValue<Number> NumberTensorValue;
typedef NumberTensorValue   Tensor;
typedef VectorValue<Number> NumberVectorValue;
typedef NumberVectorValue   Gradient;

/**
 * This class handles the computation of the L2 and/or H1 error for
 * the Systems in the EquationSystems object which is passed to it.
 *
 * \note For this to be useful, the user must attach at least one, and
 * possibly two, functions which can compute the exact solution and
 * its derivative for any component of any system.  These are the \p
 * exact_value and \p exact_deriv functions below.
 *
 * \author Benjamin S. Kirk
 * \author John W. Peterson (modifications for libmesh)
 * \date 2004
 */
class ExactSolution
{

public:
  /**
   * Constructor.  The ExactSolution object
   * must be initialized with an EquationSystems
   * object.
   */
  explicit
  ExactSolution (const EquationSystems & es);

  /**
   * The copy constructor and copy/move assignment operators are
   * deleted.  This class has containers of unique_ptrs so it can't be
   * default (shallow) copied, and it has a const reference so it
   * can't be assigned to after creation.
   */
  ExactSolution(const ExactSolution &) = delete;
  ExactSolution & operator= (const ExactSolution &) = delete;
  ExactSolution & operator= (ExactSolution &&) = delete;

  /**
   * Move constructor and destructor are defaulted out-of-line (in the
   * C file) to play nicely with our forward declarations.
   */
  ExactSolution(ExactSolution &&);
  ~ExactSolution();

  /**
   * The user can indicate that elements in certain subdomains should be
   * excluded from the error calculation by passing in a set of subdomain
   * ids to ignore. By default, all subdomains are considered in the error
   * calculation.
   */
  void set_excluded_subdomains(const std::set<subdomain_id_type> & excluded);

  /**
   * Attach function similar to system.h which
   * allows the user to attach a second EquationSystems
   * object with a reference fine grid solution.
   */
  void attach_reference_solution (const EquationSystems * es_fine);

  /**
   * Clone and attach arbitrary functors which compute the exact
   * values of the EquationSystems' solutions at any point.
   */
  void attach_exact_values (const std::vector<FunctionBase<Number> *> & f);

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
  void attach_exact_derivs (const std::vector<FunctionBase<Gradient> *> & g);

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
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.
   */
  void extra_quadrature_order (const int extraorder)
  { _extra_order = extraorder; }

  /**
   * Computes and stores the error in the solution value e = u-u_h,
   * the gradient grad(e) = grad(u) - grad(u_h), and possibly the hessian
   * grad(grad(e)) = grad(grad(u)) - grad(grad(u_h)).  Does not return
   * any value.  For that you need to call the l2_error(), h1_error()
   * or h2_error() functions respectively.
   */
  void compute_error(const std::string & sys_name,
                     const std::string & unknown_name);

  /**
   * \returns The integrated L2 error for the system \p sys_name for the
   * unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   */
  Real l2_error(const std::string & sys_name,
                const std::string & unknown_name);

  /**
   * \returns The integrated L1 error for the system \p sys_name for
   * the unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   */
  Real l1_error(const std::string & sys_name,
                const std::string & unknown_name);

  /**
   * \returns The L_INF error for the system \p sys_name for
   * the unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * compute_error() for that.
   *
   * \note The result (as for the other norms as well) is not exact,
   * but an approximation based on the chosen quadrature rule: to
   * compute it, we take the max of the absolute value of the error
   * over all the quadrature points.
   */
  Real l_inf_error(const std::string & sys_name,
                   const std::string & unknown_name);

  /**
   * \returns The H1 error for the system \p sys_name for the unknown
   * \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   */
  Real h1_error(const std::string & sys_name,
                const std::string & unknown_name);

  /**
   * \returns The H(curl) error for the system \p sys_name for the
   * unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   *
   * \note This is only valid for vector-valued elements. An error is
   * thrown if requested for scalar-valued elements.
   */
  Real hcurl_error(const std::string & sys_name,
                   const std::string & unknown_name);

  /**
   * \returns The H(div) error for the system \p sys_name for the
   * unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   *
   * \note This is only valid for vector-valued elements. An error is
   * thrown if requested for scalar-valued elements.
   */
  Real hdiv_error(const std::string & sys_name,
                  const std::string & unknown_name);

  /**
   * \returns The H2 error for the system \p sys_name for the unknown
   * \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   */
  Real h2_error(const std::string & sys_name,
                const std::string & unknown_name);

  /**
   * \returns The error in the requested norm for the system \p
   * sys_name for the unknown \p unknown_name.
   *
   * \note No error computations are actually performed, you must call
   * \p compute_error() for that.
   *
   * \note The result is not exact, but an approximation based on the
   * chosen quadrature rule.
   */
  Real error_norm(const std::string & sys_name,
                  const std::string & unknown_name,
                  const FEMNormType & norm);
private:

  /**
   * This function computes the error (in the solution and its first
   * derivative) for a single unknown in a single system.  It is a
   * private function since it is used by the implementation when
   * solving for several unknowns in several systems.
   */
  template<typename OutputShape>
  void _compute_error(const std::string & sys_name,
                      const std::string & unknown_name,
                      std::vector<Real> & error_vals);

  /**
   * This function is responsible for checking the validity of the \p
   * sys_name and \p unknown_name inputs.
   *
   * \returns A reference to the proper vector for storing the values.
   */
  std::vector<Real> & _check_inputs(const std::string & sys_name,
                                    const std::string & unknown_name);

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
   * Data structure which stores the errors:
   * ||e|| = ||u - u_h||
   * ||grad(e)|| = ||grad(u) - grad(u_h)||
   * for each unknown in a single system.
   * The name of the unknown is
   * the key for the map.
   */
  typedef std::map<std::string, std::vector<Real>> SystemErrorMap;

  /**
   * A map of SystemErrorMaps, which contains entries
   * for each system in the EquationSystems object.
   * This is required, since it is possible for two
   * systems to have unknowns with the *same name*.
   */
  std::map<std::string, SystemErrorMap> _errors;

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  const EquationSystems & _equation_systems;

  /**
   * Constant pointer to the \p EquationSystems object
   * containing the fine grid solution.
   */
  const EquationSystems * _equation_systems_fine;

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;

  /**
   * Elements in a subdomain from this set are skipped during the
   * error computation.
   */
  std::set<subdomain_id_type> _excluded_subdomains;
};



} // namespace libMesh


#endif // LIBMESH_EXACT_SOLUTION_H
