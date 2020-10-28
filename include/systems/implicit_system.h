// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_IMPLICIT_SYSTEM_H
#define LIBMESH_IMPLICIT_SYSTEM_H

// Local Includes
#include "libmesh/explicit_system.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
template <typename T> class LinearSolver;

/**
 * \brief Manages consistently variables, degrees of freedom, coefficient
 * vectors, and matrices for implicit systems.
 *
 * Implicit systems are characterized by the need to solve the (non-)linear
 * system Ax=b. This class provides, in addition to the ExplicitSystem class,
 * storage for sparse matrices. Hence this System provides means to manage a
 * right hand side vector and a sparse matrix. In addition, further matrices
 * can be managed using the ImplicitSystem::add_matrix method.
 *
 * This class provieds *no* means to solve the implicit system. This
 * functionality is provided, e.g., by the LinearImplicitSystem or the
 * NonlinearImplicitSystem class.
 *
 * \note Additional vectors/matrices can be added via parent class
 * interfaces.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class ImplicitSystem : public ExplicitSystem
{
public:

  /**
   * Constructor.
   */
  ImplicitSystem (EquationSystems & es,
                  const std::string & name,
                  const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - The destructor is defaulted out-of-line.
   */
  ImplicitSystem (const ImplicitSystem &) = delete;
  ImplicitSystem & operator= (const ImplicitSystem &) = delete;
  ImplicitSystem (ImplicitSystem &&) = default;
  ImplicitSystem & operator= (ImplicitSystem &&) = delete;
  virtual ~ImplicitSystem ();

  /**
   * The type of system.
   */
  typedef ImplicitSystem sys_type;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * The type of the parent.
   */
  typedef ExplicitSystem Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () override;

  /**
   * Prepares \p matrix and \p rhs for system assembly, then calls
   * user assembly function.
   * Can be overridden in derived classes.
   */
  virtual void assemble () override;

  /**
   * Avoids use of any cached data that might affect any solve result.  Should
   * be overridden in derived systems.
   */
  virtual void disable_cache () override;

  /**
   * \returns \p "Implicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const override { return "Implicit"; }

  /**
   * \returns A dumb pointer to a local std::unique_ptr<LinearSolver>
   * member which can be used in adjoint and/or sensitivity solves.
   *
   * To mimic the previous behavior of this function, if the
   * linear_solver member is already initialized when this function is
   * called, this function first clears it. That is, no attempt is
   * made to reuse an existing LinearSolver object. The user MUST NOT
   * attempt to clean up this pointer, otherwise the std::unique_ptr
   * held by this class will likely become corrupted.
   *
   * This function is virtual so it can be overridden in derived
   * classes if necessary, however most will probably just want to use
   * the base class behavior.
   */
  virtual LinearSolver<Number> * get_linear_solver() const;

  /**
   * \returns An integer corresponding to the upper iteration count
   * limit and a Real corresponding to the convergence tolerance to
   * be used in linear adjoint and/or sensitivity solves
   */
  virtual std::pair<unsigned int, Real>
  get_linear_solve_parameters() const;

  /**
   * Currently a no-op.
   *
   * \deprecated This function no longer needs to be called, since
   * get_linear_solver() no longer returns a heap-allocated dumb
   * pointer.
   */
  virtual void release_linear_solver(LinearSolver<Number> *) const;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   *
   * This is undefined in ImplicitSystem; subclasses each have their
   * own way of handling assembly.
   */
  virtual void assembly (bool /* get_residual */,
                         bool /* get_jacobian */,
                         bool /* apply_heterogeneous_constraints */ = false,
                         bool /* apply_no_constraints */ = false)
  { libmesh_not_implemented(); }

  /**
   * Residual parameter derivative function.
   *
   * Uses finite differences by default.
   *
   * This will assemble the sensitivity rhs vectors to hold
   * -(partial R / partial p_i), making them ready to solve
   * the forward sensitivity equation.
   *
   * Can be overridden in derived classes.
   */
  virtual void assemble_residual_derivatives (const ParameterVector & parameters) override;

  /**
   * For explicit systems, just assemble and solve the system A*x=b.
   * Should be overridden in derrived systems to provide a solver for the
   * system.
   */
  virtual void solve () override
  { libmesh_not_implemented(); }

  /**
   * Assembles & solves the linear system(s) (dR/du)*u_p = -dR/dp, for
   * those parameters contained within \p parameters.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
  sensitivity_solve (const ParameterVector & parameters) override;

  /**
   * Assembles & solves the linear system(s) (dR/du)*u_w = sum(w_p*-dR/dp), for
   * those parameters p contained within \p parameters weighted by the
   * values w_p found within \p weights.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
  weighted_sensitivity_solve (const ParameterVector & parameters,
                              const ParameterVector & weights) override;

  /**
   * Assembles & solves the linear system (dR/du)^T*z = dq/du, for
   * those quantities of interest q specified by \p qoi_indices.
   *
   * Leave \p qoi_indices empty to solve all adjoint problems.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
  adjoint_solve (const QoISet & qoi_indices = QoISet()) override;

  /**
   * Assembles & solves the linear system(s)
   * (dR/du)^T*z_w = sum(w_p*(d^2q/dudp - d^2R/dudp*z)), for those
   * parameters p contained within \p parameters, weighted by the
   * values w_p found within \p weights.
   *
   * Assumes that adjoint_solve has already calculated z for each qoi
   * in \p qoi_indices.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
  weighted_sensitivity_adjoint_solve (const ParameterVector & parameters,
                                      const ParameterVector & weights,
                                      const QoISet & qoi_indices = QoISet()) override;

  /**
   * Solves for the derivative of each of the system's quantities of
   * interest q in \p qoi[qoi_indices] with respect to each parameter in
   * \p parameters, placing the result for qoi \p i and parameter \p j
   * into \p sensitivities[i][j].
   *
   * Uses adjoint_solve() and the adjoint sensitivity method.
   *
   * Currently uses finite differenced derivatives (partial q /
   * partial p) and (partial R / partial p).
   */
  virtual void adjoint_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                  const ParameterVector & parameters,
                                                  SensitivityData & sensitivities) override;

  /**
   * Solves for the derivative of each of the system's quantities of
   * interest q in \p qoi[qoi_indices] with respect to each parameter in
   * \p parameters, placing the result for qoi \p i and parameter \p j
   * into \p sensitivities[i][j].
   *
   * Uses the forward sensitivity method.
   *
   * Currently uses finite differenced derivatives (partial q /
   * partial p) and (partial R / partial p).
   */
  virtual void forward_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                  const ParameterVector & parameters,
                                                  SensitivityData & sensitivities) override;

  /**
   * For each of the system's quantities of interest q in
   * \p qoi[qoi_indices], and for a vector of parameters p, the
   * parameter sensitivity Hessian H_ij is defined as
   * H_ij = (d^2 q)/(d p_i d p_j)
   * This Hessian is the output of this method, where for each q_i,
   * H_jk is stored in \p hessian.second_derivative(i,j,k).
   *
   * Note that in some cases only
   * \link current_local_solution \endlink is used during assembly,
   * and, therefore, if \link solution \endlink has been altered
   * without \link update() \endlink being called, then the
   * user must call \link update() \endlink before calling
   * this function.
   */
  virtual void qoi_parameter_hessian(const QoISet & qoi_indices,
                                     const ParameterVector & parameters,
                                     SensitivityData & hessian) override;

  /**
   * For each of the system's quantities of interest q in
   * \p qoi[qoi_indices], and for a vector of parameters p, the
   * parameter sensitivity Hessian H_ij is defined as
   * H_ij = (d^2 q)/(d p_i d p_j)
   * The Hessian-vector product, for a vector v_k in parameter space, is
   * S_j = H_jk v_k
   * This product is the output of this method, where for each q_i,
   * S_j is stored in \p sensitivities[i][j].
   */
  virtual void qoi_parameter_hessian_vector_product(const QoISet & qoi_indices,
                                                    const ParameterVector & parameters,
                                                    const ParameterVector & vector,
                                                    SensitivityData & product) override;

  /**
   * \returns A const reference to the system's primary matrix.
   */
  const SparseMatrix<Number> & get_system_matrix() const;

  /**
   * \returns A reference to the system's primary matrix.
   */
  SparseMatrix<Number> & get_system_matrix();

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * system matrix A.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_system_matrix() instead.
   */
  SparseMatrix<Number> * matrix;

  /**
   * By default, the system will zero out the matrix and the right hand side.
   * If this flag is false, it is the responsibility of the client code
   * to take care of setting these to zero before assembly begins
   */
  bool zero_out_matrix_and_rhs;

  /**
   * This class handles all the details of interfacing with various
   * linear algebra packages like PETSc or LASPACK.  This is a public
   * member for backwards compatibility reasons, but in general it's
   * better to use get_linear_solver() to access this member, since
   * that function will also handle initialization if it hasn't
   * already been taken care of.
   */
  mutable std::unique_ptr<LinearSolver<Number>> linear_solver;

protected:
  /**
   * Adds the system matrix
   */
  virtual void add_matrices() override;
};


} // namespace libMesh

#endif // LIBMESH_IMPLICIT_SYSTEM_H
