// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TRILINOS_AZTEC_LINEAR_SOLVER_H
#define LIBMESH_TRILINOS_AZTEC_LINEAR_SOLVER_H



// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/linear_solver.h"

/**
 * Trilinos include files.
 */
#ifdef LIBMESH_TRILINOS_HAVE_AZTECOO
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * This class provides an interface to AztecOO
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * \author Benjamin Kirk
 * \date 2002-2008
 */
template <typename T>
class AztecLinearSolver : public LinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Aztec data structures
   */
  AztecLinearSolver (const libMesh::Parallel::Communicator & comm
                     LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  ~AztecLinearSolver ();

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () libmesh_override;

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init (const char * name=libmesh_nullptr) libmesh_override;

  /**
   * Call the Aztec solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real>
  solve (SparseMatrix<T> & matrix_in,
         NumericVector<T> & solution_in,
         NumericVector<T> & rhs_in,
         const double tol,
         const unsigned int m_its) libmesh_override
  {
    return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
  }

  /**
   * This method allows you to call a linear solver while specifying
   * the matrix to use as the (left) preconditioning matrix.  Note
   * that the linear solver will not compute a preconditioner in this
   * case, and will instead premultiply by the matrix you provide.
   */
  virtual std::pair<unsigned int, Real>
  solve (SparseMatrix<T> & matrix,
         SparseMatrix<T> & preconditioner,
         NumericVector<T> & solution,
         NumericVector<T> & rhs,
         const double tol,
         const unsigned int m_its) libmesh_override;

  /**
   * This function solves a system whose matrix is a shell matrix.
   */
  virtual std::pair<unsigned int, Real>
  solve (const ShellMatrix<T> & shell_matrix,
         NumericVector<T> & solution_in,
         NumericVector<T> & rhs_in,
         const double tol,
         const unsigned int m_its) libmesh_override;

  /**
   * This function solves a system whose matrix is a shell matrix, but
   * a sparse matrix is used as preconditioning matrix, this allowing
   * other preconditioners than JACOBI.
   */
  virtual std::pair<unsigned int, Real>
  solve (const ShellMatrix<T> & shell_matrix,
         const SparseMatrix<T> & precond_matrix,
         NumericVector<T> & solution_in,
         NumericVector<T> & rhs_in,
         const double tol,
         const unsigned int m_its) libmesh_override;

  /**
   * Fills the input vector with the sequence of residual norms
   * from the latest iterative solve.
   */
  void get_residual_history(std::vector<double> & hist);

  /**
   * Returns just the initial residual for the solve just
   * completed with this interface.  Use this method instead
   * of the one above if you just want the starting residual
   * and not the entire history.
   */
  Real get_initial_residual();

  /**
   * Returns the solver's convergence flag
   */
  virtual LinearConvergenceReason get_converged_reason() const libmesh_override;

private:

  /**
   * Tells AztecOO to use the user-specified solver stored in
   * \p _solver_type
   */
  void set_solver_type ();

  /**
   * The Epetra linear problem object.
   */
  Epetra_LinearProblem * _linear_problem;

  /**
   * The AztecOO solver object
   */
  AztecOO * _linear_solver;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
AztecLinearSolver<T>::AztecLinearSolver (const libMesh::Parallel::Communicator & comm) :
  LinearSolver<T>(comm)
{
  if (this->n_processors() == 1)
    this->_preconditioner_type = ILU_PRECOND;
  else
    this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
}



template <typename T>
inline
AztecLinearSolver<T>::~AztecLinearSolver ()
{
  this->clear ();
}

} // namespace libMesh



#endif // #ifdef LIBMESH_TRILINOS_HAVE_AZTECOO
#endif // LIBMESH_TRILINOS_AZTEC_LINEAR_SOLVER_H
