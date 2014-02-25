// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EIGEN_SPARSE_LINEAR_SOLVER_H
#define LIBMESH_EIGEN_SPARSE_LINEAR_SOLVER_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// Eigen includes

// Local includes
#include "libmesh/linear_solver.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/eigen_sparse_matrix.h"

// C++ includes

namespace libMesh
{



/**
 * This class provides an interface to Eigen
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * @author Benjamin Kirk, 2013
 */
template <typename T>
class EigenSparseLinearSolver : public LinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Eigen data structures
   */
  EigenSparseLinearSolver (const libMesh::Parallel::Communicator &comm
                           LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  ~EigenSparseLinearSolver ();

  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();

  /**
   * Call the Eigen solver
   */
  std::pair<unsigned int, Real>
  solve (SparseMatrix<T>  &matrix,
         NumericVector<T> &solution,
         NumericVector<T> &rhs,
         const double tol,
         const unsigned int m_its);

  /**
   * Call the Eigen solver to solve A^T x = b
   */
  std::pair<unsigned int, Real>
  adjoint_solve (SparseMatrix<T>  &matrix,
                 NumericVector<T> &solution,
                 NumericVector<T> &rhs,
                 const double tol,
                 const unsigned int m_its);

  /**
   * Call the Eigen solver
   */
  std::pair<unsigned int, Real>
  solve (SparseMatrix<T>  &matrix,
         SparseMatrix<T>  &pc,
         NumericVector<T> &solution,
         NumericVector<T> &rhs,
         const double tol,
         const unsigned int m_its);

  /**
   * This function solves a system whose matrix is a shell matrix.
   */
  std::pair<unsigned int, Real>
  solve (const ShellMatrix<T>& shell_matrix,
         NumericVector<T>& solution_in,
         NumericVector<T>& rhs_in,
         const double tol,
         const unsigned int m_its);

  /**
   * This function solves a system whose matrix is a shell matrix, but
   * a sparse matrix is used as preconditioning matrix, this allowing
   * other preconditioners than JACOBI.
   */
  virtual std::pair<unsigned int, Real>
  solve (const ShellMatrix<T>& shell_matrix,
         const SparseMatrix<T>& precond_matrix,
         NumericVector<T>& solution_in,
         NumericVector<T>& rhs_in,
         const double tol,
         const unsigned int m_its);

  /**
   * Prints a useful message about why the latest linear solve
   * con(di)verged.
   */
  virtual void print_converged_reason() const;

  /**
   * Returns the solver's convergence flag
   */
  virtual LinearConvergenceReason get_converged_reason() const;

private:

  /**
   * Tells Eigen to use the user-specified preconditioner stored in
   * \p _preconditioner_type
   */
  void set_eigen_preconditioner_type ();

};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
EigenSparseLinearSolver<T>::EigenSparseLinearSolver (const libMesh::Parallel::Communicator &comm) :
  LinearSolver<T>(comm)
{
}



template <typename T>
inline
EigenSparseLinearSolver<T>::~EigenSparseLinearSolver ()
{
  this->clear ();
}



template <typename T>
inline
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (SparseMatrix<T>&,
                                   SparseMatrix<T>&,
                                   NumericVector<T>&,
                                   NumericVector<T>&,
                                   const double,
                                   const unsigned int)
{
  libMesh::err << "ERROR: Eigen does not support a user-supplied preconditioner!"
               << std::endl;
  libmesh_error();

  std::pair<unsigned int, Real> p;
  return p;
}

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_SPARSE_LINEAR_SOLVER_H
