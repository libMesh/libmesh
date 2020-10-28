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


#ifndef LIBMESH_EIGEN_SYSTEM_H
#define LIBMESH_EIGEN_SYSTEM_H

#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// Local Includes
#include "libmesh/system.h"
#include "libmesh/eigen_solver.h"

namespace libMesh
{

// Forward Declarations
template <typename T> class SparseMatrix;
template <typename T> class ShellMatrix;


/**
 * \brief Manages consistently variables, degrees of freedom, and coefficient
 * vectors for eigenvalue problems.
 *
 * Currently, this class is able to handle standard eigenvalue
 * problems \p A*x=lambda*x and generalized eigenvalue problems \p
 * A*x=lambda*B*x.
 *
 * The matrices EigenSystem::matrix_A and EigenSystem::matrix_B should be
 * filled during assembly.
 *
 * \author Steffen Peterson
 * \date 2005
 * \brief Base class for defining systems of equations for eigenproblems.
 */
class EigenSystem : public System
{
public:

  /**
   * Constructor.
   */
  EigenSystem (EquationSystems & es,
               const std::string & name_in,
               const unsigned int number_in);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - Destructor is defaulted out-of-line
   */
  EigenSystem (const EigenSystem &) = delete;
  EigenSystem & operator= (const EigenSystem &) = delete;
  EigenSystem (EigenSystem &&) = default;
  EigenSystem & operator= (EigenSystem &&) = delete;
  virtual ~EigenSystem ();

  /**
   * The type of system.
   */
  typedef EigenSystem sys_type;

  /**
   * The type of the parent
   */
  typedef System Parent;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

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
   * Assembles & solves the eigen system.
   */
  virtual void solve () override;

  /**
   * \returns Real and imaginary part of the ith eigenvalue and copies
   * the respective eigen vector to the solution vector.
   */
  virtual std::pair<Real, Real> get_eigenpair (dof_id_type i);

  /**
   * \returns Real and imaginary part of the ith eigenvalue but
   * does not copy the respective eigen vector to the solution vector.
   */
  virtual std::pair<Real, Real> get_eigenvalue (dof_id_type i);

  /**
   * \returns \p "Eigen".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const override { return "Eigen"; }

  /**
   * \returns The number of converged eigenpairs.
   */
  unsigned int get_n_converged () const {return _n_converged_eigenpairs;}

  /**
   * \returns The number of eigen solver iterations.
   */
  unsigned int get_n_iterations () const {return _n_iterations;}

  /**
   * Sets the type of the current eigen problem.
   */
  void set_eigenproblem_type (EigenProblemType ept);

  /**
   * \returns The eigen problem type.
   */
  EigenProblemType get_eigenproblem_type () const {return _eigen_problem_type;}

  /**
   * Sets an initial eigen vector
   */
  void set_initial_space(NumericVector<Number> & initial_space_in);

  /**
   * \returns \p true if the underlying problem is generalized
   * , false otherwise.
   */
  bool generalized () const;

  /**
   * \returns \p true if the shell matrices are used
   */
  bool use_shell_matrices() const { return _use_shell_matrices; }

  /**
   * Set a flag to use shell matrices
   */
  void use_shell_matrices(bool use_shell_matrices) { _use_shell_matrices = use_shell_matrices; }

  /**
   * \returns \p true if a shell preconditioning matrix is used
   */
  bool use_shell_precond_matrix() const { return _use_shell_precond_matrix; }

  /**
   * Set a flag to use a shell preconditioning matrix
   */
  void use_shell_precond_matrix(bool use_shell_precond_matrix) { _use_shell_precond_matrix = use_shell_precond_matrix; }

  /**
   * \returns A const reference to the system matrix used for standard eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices.
   */
  const SparseMatrix<Number> & get_matrix_A() const;

  /**
   * \returns A reference to the system matrix used for standard eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices.
   */
  SparseMatrix<Number> & get_matrix_A();

  /**
   * \returns A const reference to the system matrix used for generalized eigenvalue problems
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices and generalized().
   */
  const SparseMatrix<Number> & get_matrix_B() const;

  /**
   * \returns A const reference to the system matrix used for generalized eigenvalue problems
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices and generalized().
   */
  SparseMatrix<Number> & get_matrix_B();

  /**
   * \returns A const reference to the preconditioning matrix.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices.
   */
  const SparseMatrix<Number> & get_precond_matrix() const;

  /**
   * \returns A reference to the preconditioning matrix.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if !_use_shell_matrices.
   */
  SparseMatrix<Number> & get_precond_matrix();

  /**
   * \returns A const reference to the system shell matrix used for standard eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if _use_shell_matrices.
   */
  const ShellMatrix<Number> & get_shell_matrix_A() const;

  /**
   * \returns A reference to the system shell matrix used for standard eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if _use_shell_matrices.
   */
  ShellMatrix<Number> & get_shell_matrix_A();

  /**
   * \returns A const reference to the system shell matrix used for generalized eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if _use_shell_matrices and generalized().
   */
  const ShellMatrix<Number> & get_shell_matrix_B() const;

  /**
   * \returns A reference to the system shell matrix used for generalized eigenvalue problems.
   *
   * This matrix should only be available (and therefore this should only be called)
   * if _use_shell_matrices and generalized().
   */
  ShellMatrix<Number> & get_shell_matrix_B();

  /**
   * \returns A const reference to the system shell matrix used for preconditioning.
   */
  const ShellMatrix<Number> & get_shell_precond_matrix() const;

  /**
   * \returns A reference to the system shell matrix used for preconditioning.
   */
  ShellMatrix<Number> & get_shell_precond_matrix();

  /**
   * \returns A const reference to the EigenSolver.
   */
  const EigenSolver<Number> & get_eigen_solver() const;
  /**
   * \returns A reference to the EigenSolver.
   */
  EigenSolver<Number> & get_eigen_solver();

  /**
   * \returns Whether or not the system has matrix A
   */
  bool has_matrix_A() const;

  /**
   * \returns Whether or not the system has matrix B
   */
  bool has_matrix_B() const;

  /**
   * \returns Whether or not the system has the non-shell preconditioning matrix
   */
  bool has_precond_matrix() const;

  /**
   * \returns Whether or not the system has the shell matrix A
   */
  bool has_shell_matrix_A() const;

  /**
   * \returns Whether or not the system has the shell matrix B
   */
  bool has_shell_matrix_B() const;

  /**
   * \returns Whether or not the system has the shell preconditioning matrix
   */
  bool has_shell_precond_matrix() const;

  /**
   * The system matrix for standard eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_matrix_A() instead.
   */
  SparseMatrix<Number> * matrix_A;

  /**
   * A second system matrix for generalized eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_matrix_B() instead.
   */
  SparseMatrix<Number> * matrix_B;

  /**
   * The system shell matrix for standard eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_shell_matrix_A() instead.
   */
  std::unique_ptr<ShellMatrix<Number>> shell_matrix_A;

  /**
   * A second system shell matrix for generalized eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_shell_matrix_B() instead.
   */
  std::unique_ptr<ShellMatrix<Number>> shell_matrix_B;

  /**
   * A preconditioning matrix
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_precond_matrix() instead.
   */
  SparseMatrix<Number> * precond_matrix;

  /**
   * A preconditioning shell matrix
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_shell_precond_matrix() instead.
   */
  std::unique_ptr<ShellMatrix<Number>> shell_precond_matrix;

  /**
   * The EigenSolver, defining which interface, i.e solver
   * package to use.
   *
   * Public access to this member variable will be depreacted in
   * the future! Use get_eigen_solver instead.
   */
  std::unique_ptr<EigenSolver<Number>> eigen_solver;


protected:

  /**
   * Adds the necessary matrices and shell matrices
   */
  virtual void add_matrices () override;

  /**
   * Initializes the matrices associated with the system
   */
  virtual void init_matrices () override;

  /**
   * Set the _n_converged_eigenpairs member, useful for
   * subclasses of EigenSystem.
   */
  void set_n_converged (unsigned int nconv)
  { _n_converged_eigenpairs = nconv; }

  /**
   * Set the _n_iterations member, useful for subclasses of
   * EigenSystem.
   */
  void set_n_iterations (unsigned int its)
  { _n_iterations = its;}


private:

  /**
   * The number of converged eigenpairs.
   */
  unsigned int _n_converged_eigenpairs;

  /**
   * The number of iterations of the eigen solver algorithm.
   */
  unsigned int _n_iterations;

  /**
   * The type of the eigenvalue problem.
   */
  EigenProblemType _eigen_problem_type;

  /**
   * A boolean flag to indicate whether or not to use shell matrices
   */
  bool _use_shell_matrices;

  /**
   * A boolean flag to indicate whether or not to use a shell preconditioning matrix
   */
  bool _use_shell_precond_matrix;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC

#endif // LIBMESH_EIGEN_SYSTEM_H
