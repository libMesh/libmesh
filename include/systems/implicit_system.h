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
template <typename T> class SparseMatrix;

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
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ImplicitSystem (EquationSystems & es,
                  const std::string & name,
                  const unsigned int number);

  /**
   * Destructor.
   */
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
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () override;

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
   * \returns A pointer to a linear solver appropriate for use in
   * adjoint and/or sensitivity solves
   *
   * This function must be overridden in derived classes, since this
   * base class does not have a valid LinearSolver to hand back a
   * pointer to.
   *
   * \deprecated This function's current behavior, i.e. allocating a
   * LinearSolver and handing it back to the user, makes it very easy
   * to leak memory, and probably won't have the intended effect,
   * i.e. of setting some parameters on a LinearSolver that the System
   * would later use internally.
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
   * Releases a pointer to a linear solver acquired by
   * \p this->get_linear_solver()
   *
   * \deprecated This function is designed to work with the deprecated
   * get_linear_solver() function, so its use is now deprecated as
   * well.
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
   * Matrix iterator typedefs.
   */
  typedef std::map<std::string, SparseMatrix<Number> *>::iterator        matrices_iterator;
  typedef std::map<std::string, SparseMatrix<Number> *>::const_iterator  const_matrices_iterator;

  /**
   * Adds the additional matrix \p mat_name to this system.  Only
   * allowed prior to \p assemble().  All additional matrices
   * have the same sparsity pattern as the matrix used during
   * solution.  When not \p System but the user wants to
   * initialize the mayor matrix, then all the additional matrices,
   * if existent, have to be initialized by the user, too.
   *
   * This non-template method will add a derived matrix type corresponding to
   * the solver package. If the user wishes to specify the matrix type to add,
   * use the templated \p add_matrix method instead
   *
   * @param type The serial/parallel/ghosted type of the matrix
   *
   */
  SparseMatrix<Number> & add_matrix (const std::string & mat_name, ParallelType type = PARALLEL);

  /**
   * Adds the additional matrix \p mat_name to this system.  Only
   * allowed prior to \p assemble().  All additional matrices
   * have the same sparsity pattern as the matrix used during
   * solution.  When not \p System but the user wants to
   * initialize the mayor matrix, then all the additional matrices,
   * if existent, have to be initialized by the user, too.
   *
   * This method will create add a derived matrix of type
   * \p MatrixType<Number>. One can use the non-templated \p add_matrix method to
   * add a matrix corresponding to the default solver package
   *
   * @param type The serial/parallel/ghosted type of the matrix
   */
  template <template <typename> class>
  SparseMatrix<Number> & add_matrix (const std::string & mat_name, ParallelType = PARALLEL);

  /**
   * Removes the additional matrix \p mat_name from this system
   */
  void remove_matrix(const std::string & mat_name);

  /**
   * \returns \p true if this \p System has a matrix associated with the
   * given name, \p false otherwise.
   */
  bool have_matrix (const std::string & mat_name) const;

  /**
   * \returns A const pointer to this system's additional matrix
   * named \p mat_name, or \p nullptr if no matrix by that name
   * exists.
   */
  const SparseMatrix<Number> * request_matrix (const std::string & mat_name) const;

  /**
   * \returns A writable pointer to this system's additional matrix
   * named \p mat_name, or \p nullptr if no matrix by that name
   * exists.
   */
  SparseMatrix<Number> * request_matrix (const std::string & mat_name);

  /**
   * \returns A const reference to this system's additional matrix
   * named \p mat_name.
   *
   * None of these matrices is involved in the solution process.
   * Access is only granted when the matrix is already properly
   * initialized.
   */
  const SparseMatrix<Number> & get_matrix (const std::string & mat_name) const;

  /**
   * \returns A writable reference to this system's additional matrix
   * named \p mat_name.
   *
   * None of these matrices is involved in the solution process.
   * Access is only granted when the matrix is already properly
   * initialized.
   */
  SparseMatrix<Number> & get_matrix (const std::string & mat_name);

  /**
   * \returns The number of matrices handled by this system
   */
  virtual unsigned int n_matrices () const override;

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * system matrix A.
   */
  SparseMatrix<Number> * matrix;

  /**
   * By default, the system will zero out the matrix and the right hand side.
   * If this flag is false, it is the responsibility of the client code
   * to take care of setting these to zero before assembly begins
   */
  bool zero_out_matrix_and_rhs;

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data () override;

  /**
   * Initializes the matrices associated with this system.
   */
  virtual void init_matrices ();

private:

  /**
   * Add the system matrix to the \p _matrices data structure.
   * Useful in initialization.
   */
  void add_system_matrix ();

  /**
   * Some systems need an arbitrary number of matrices.
   */
  std::map<std::string, SparseMatrix<Number> *> _matrices;

  /**
   * Holds the types of the matrices
   */
  std::map<std::string, ParallelType> _matrix_types;

  /**
   * \p true when additional matrices may still be added, \p false otherwise.
   */
  bool _can_add_matrices;
};



// ------------------------------------------------------------
// ImplicitSystem inline methods
inline
bool ImplicitSystem::have_matrix (const std::string & mat_name) const
{
  return (_matrices.count(mat_name));
}


inline
unsigned int ImplicitSystem::n_matrices () const
{
  return cast_int<unsigned int>(_matrices.size());
}

template <template <typename> class MatrixType>
inline
SparseMatrix<Number> &
ImplicitSystem::add_matrix (const std::string & mat_name,
                            const ParallelType type)
{
  // only add matrices before initializing...
  libmesh_error_msg_if(!_can_add_matrices,
                       "ERROR: Too late.  Cannot add matrices to the system after initialization"
                       "\n any more.  You should have done this earlier.");

  // Return the matrix if it is already there.
  if (this->have_matrix(mat_name))
    return *(_matrices[mat_name]);

  // Otherwise build the matrix and return it.
  SparseMatrix<Number> * buf = libmesh_make_unique<MatrixType<Number>>(this->comm()).release();
  _matrices.emplace(mat_name, buf);
  _matrix_types.emplace(mat_name, type);

  return *buf;
}

} // namespace libMesh

#endif // LIBMESH_IMPLICIT_SYSTEM_H
