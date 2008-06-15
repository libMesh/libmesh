// $Id: trilinos_azec_linear_solver.h 2501 2007-11-20 02:33:29Z benkirk $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __trilinos_azec_linear_solver_h__
#define __trilinos_azec_linear_solver_h__



// Local includes
#include "linear_solver.h"

/**
 * Trilinos include files.
 */
#ifdef HAVE_TRILINOS
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

// C++ includes
#include <vector>

/**
 * This class provides an interface to AztecOO
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2008
 */

template <typename T>
class AztecLinearSolver : public LinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Petsc data structures
   */
  AztecLinearSolver ();
    
  /**
   * Destructor.
   */
  ~AztecLinearSolver ();
  
  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();

  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */    
  std::pair<unsigned int, Real> 
  solve (SparseMatrix<T>  &matrix_in,
	 NumericVector<T> &solution_in,
	 NumericVector<T> &rhs_in,
	 const double tol,
	 const unsigned int m_its)
  {
    return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
  }

  /**
   * This method allows you to call a linear solver while specifying
   * the matrix to use as the (left) preconditioning matrix.  Note
   * that the linear solver will not compute a preconditioner in this
   * case, and will instead premultiply by the matrix you provide.
   *
   * In PETSc, this is accomplished by calling
   *
   * PCSetType(_pc, PCMAT);
   *
   * before invoking KSPSolve().  Note: this functionality is not implemented
   * in the LinearSolver class since there is not a built-in analog
   * to this method for LasPack -- You could probably implement it by hand
   * if you wanted.
   */
  std::pair<unsigned int, Real> 
  solve (SparseMatrix<T>  &matrix,
	 SparseMatrix<T>  &preconditioner,
	 NumericVector<T> &solution,
	 NumericVector<T> &rhs,
	 const double tol,
	 const unsigned int m_its);
  
  /**
   * Fills the input vector with the sequence of residual norms
   * from the latest iterative solve.
   */
  void get_residual_history(std::vector<double>& hist);

  /**
   * Returns just the initial residual for the solve just
   * completed with this interface.  Use this method instead
   * of the one above if you just want the starting residual
   * and not the entire history.
   */
  Real get_initial_residual();

  /**
   * Prints a useful message about why the latest linear solve
   * con(di)verged.
   */
  virtual void print_converged_reason();
  
private:

  /**
   * The Epetra linear problem object.
   */
  Epetra_LinearProblem _linear_problem;

  /**
   * The AztecOO solver object
   */
  AztecOO _linear_solver;
  

};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
AztecLinearSolver<T>::AztecLinearSolver ()
{
  LIBMESH_THROW(libMesh::NotImplemented());

  if (libMesh::n_processors() == 1)
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



#endif // #ifdef HAVE_TRILINOS
#endif // #ifdef __trilinos_azec_linear_solver_h__
