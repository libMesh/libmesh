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



#ifndef LIBMESH_IMPLICIT_SYSTEM_H
#define LIBMESH_IMPLICIT_SYSTEM_H

// Local Includes
#include "libmesh/explicit_system.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
template <typename T> class LinearSolver;
template <typename T> class SparseMatrix;



/**
 * This class provides a specific system class.  It aims
 * at implicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p ExplicitSystem.
 */

// ------------------------------------------------------------
// ImplicitSystem class definition

class ImplicitSystem : public ExplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ImplicitSystem (EquationSystems& es,
		  const std::string& name,
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
   * @returns a clever pointer to the system.
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
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();

  /**
   * Prepares \p matrix and \p rhs for system assembly, then calls
   * user assembly function.
   * @e Can be overloaded in derived classes.
   */
  virtual void assemble ();

//   /**
//    * Assembles & solves the linear system Ax=b.
//    */
//   virtual void solve ();


  /**
   * Avoids use of any cached data that might affect any solve result.  Should
   * be overloaded in derived systems.
   */
  virtual void disable_cache ();

  /**
   * @returns \p "Implicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Implicit"; }

  /**
   * Returns a pointer to a linear solver appropriate for use in
   * adjoint and/or sensitivity solves
   *
   * In \p ImplicitSystem this just gives us a default LinearSolver;
   * but it may be overridden in derived systems, e.g. to set solver
   * parameters
   */
  virtual LinearSolver<Number> *get_linear_solver() const;

  /**
   * Returns an integer corresponding to the upper iteration count
   * limit and a Real corresponding to the convergence tolerance to
   * be used in linear adjoint and/or sensitivity solves
   */
  virtual std::pair<unsigned int, Real>
    get_linear_solve_parameters() const;

  /**
   * Releases a pointer to a linear solver acquired by
   * \p this->get_linear_solver()
   */
  virtual void release_linear_solver(LinearSolver<Number> *) const;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   *
   * This is undefined in ImplicitSystem; subclasses each have their
   * own way of handling assembly.
   */
  virtual void assembly(bool /* get_residual */ , bool /* get_jacobian */)
    { libmesh_error(); }

  /**
   * Residual parameter derivative function.
   *
   * Uses finite differences by default.
   *
   * This will assemble the sensitivity rhs vectors to hold
   * -(partial R / partial p_i), making them ready to solve
   * the forward sensitivity equation.
   *
   * @e Can be overloaded in derived classes.
   */
  virtual void assemble_residual_derivatives (const ParameterVector& parameters);

  /**
   * Assembles & solves the linear system(s) (dR/du)*u_p = -dR/dp, for
   * those parameters contained within \p parameters.
   *
   * Returns a pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
    sensitivity_solve (const ParameterVector& parameters);

  /**
   * Assembles & solves the linear system(s) (dR/du)*u_w = sum(w_p*-dR/dp), for
   * those parameters p contained within \p parameters weighted by the
   * values w_p found within \p weights.
   *
   * Returns a pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
    weighted_sensitivity_solve (const ParameterVector& parameters,
                                const ParameterVector& weights);

  /**
   * Assembles & solves the linear system (dR/du)^T*z = dq/du, for
   * those quantities of interest q specified by \p qoi_indices.
   *
   * Leave \p qoi_indices empty to solve all adjoint problems.
   *
   * Returns a pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
    adjoint_solve (const QoISet& qoi_indices = QoISet());

  /**
   * Assembles & solves the linear system(s)
   * (dR/du)^T*z_w = sum(w_p*(d^2q/dudp - d^2R/dudp*z)), for those
   * parameters p contained within \p parameters, weighted by the
   * values w_p found within \p weights.
   *
   * Assumes that adjoint_solve has already calculated z for each qoi
   * in \p qoi_indices.
   *
   * Returns a pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   */
  virtual std::pair<unsigned int, Real>
    weighted_sensitivity_adjoint_solve (const ParameterVector& parameters,
                                        const ParameterVector& weights,
                                        const QoISet& qoi_indices = QoISet());

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
  virtual void adjoint_qoi_parameter_sensitivity (const QoISet& qoi_indices,
                                                  const ParameterVector& parameters,
                                                  SensitivityData& sensitivities);

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
  virtual void forward_qoi_parameter_sensitivity (const QoISet& qoi_indices,
                                                  const ParameterVector& parameters,
                                                  SensitivityData& sensitivities);

  /**
   * For each of the system's quantities of interest q in
   * \p qoi[qoi_indices], and for a vector of parameters p, the
   * parameter sensitivity Hessian H_ij is defined as
   * H_ij = (d^2 q)/(d p_i d p_j)
   * This Hessian is the output of this method, where for each q_i,
   * H_jk is stored in \p hessian.second_derivative(i,j,k).
   */
  virtual void qoi_parameter_hessian(const QoISet& qoi_indices,
                                     const ParameterVector& parameters,
                                     SensitivityData& hessian);

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
  virtual void qoi_parameter_hessian_vector_product(const QoISet& qoi_indices,
                                                    const ParameterVector& parameters,
                                                    const ParameterVector& vector,
                                                    SensitivityData& product);

  /**
   * Matrix iterator typedefs.
   */
  typedef std::map<std::string, SparseMatrix<Number>* >::iterator        matrices_iterator;
  typedef std::map<std::string, SparseMatrix<Number>* >::const_iterator  const_matrices_iterator;

  /**
   * Adds the additional matrix \p mat_name to this system.  Only
   * allowed @e prior to \p assemble().  All additional matrices
   * have the same sparsity pattern as the matrix used during
   * solution.  When not \p System but the @e user wants to
   * initialize the mayor matrix, then all the additional matrices,
   * if existent, have to be initialized by the user, too.
   */
  SparseMatrix<Number> & add_matrix (const std::string& mat_name);

  /**
   * @returns \p true if this \p System has a matrix associated with the
   * given name, \p false otherwise.
   */
  bool have_matrix (const std::string& mat_name) const;

  /**
   * @returns a const pointer to this system's @e additional matrix
   * named \p mat_name, or returns \p NULL if no matrix by that name
   * exists.
   */
  const SparseMatrix<Number> * request_matrix (const std::string& mat_name) const;

  /**
   * @returns a writable pointer to this system's @e additional matrix
   * named \p mat_name, or returns \p NULL if no matrix by that name
   * exists.
   */
  SparseMatrix<Number> * request_matrix (const std::string& mat_name);

  /**
   * @returns a const reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  const SparseMatrix<Number> & get_matrix (const std::string& mat_name) const;

  /**
   * @returns a writeable reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  SparseMatrix<Number> & get_matrix (const std::string& mat_name);

  /**
   * @returns the number of matrices handled by this system
   */
  virtual unsigned int n_matrices () const;

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * system matrix A.
   */
  SparseMatrix<Number> * matrix;



protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

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
  std::map<std::string, SparseMatrix<Number>* > _matrices;

  /**
   * \p true when additional matrices may still be added, \p false otherwise.
   */
  bool _can_add_matrices;
};



// ------------------------------------------------------------
// ImplicitSystem inline methods
inline
bool ImplicitSystem::have_matrix (const std::string& mat_name) const
{
  return (_matrices.count(mat_name));
}


inline
unsigned int ImplicitSystem::n_matrices () const
{
 return _matrices.size();
}

} // namespace libMesh

#endif // LIBMESH_IMPLICIT_SYSTEM_H
