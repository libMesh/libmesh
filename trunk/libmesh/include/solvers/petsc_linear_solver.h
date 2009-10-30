// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __petsc_linear_solver_h__
#define __petsc_linear_solver_h__

#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// C++ includes

// Local includes
#include "linear_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_macro.h"

/**
 * Petsc include files.
 */

EXTERN_C_FOR_PETSC_BEGIN
#if PETSC_VERSION_LESS_THAN(2,2,0)
#  include <petscsles.h>
#else
#  include <petscksp.h>
#endif
EXTERN_C_FOR_PETSC_END

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed for preconditioning
// 
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  /**
   * This function is called by PETSc to initialize the preconditioner.
   * ctx will hold the Preconditioner.
   */
  PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx);

  /**
   * This function is called by PETSc to acctually apply the preconditioner.
   * ctx will hold the Preconditioner.
   */
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y);
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC);
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC, Vec x, Vec y);
#endif
} // end extern "C"


/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2007
 */

template <typename T>
class PetscLinearSolver : public LinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Petsc data structures
   */
  PetscLinearSolver ();
    
  /**
   * Destructor.
   */
  ~PetscLinearSolver ();
  
  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();
  
  /**
   * Initialize data structures if not done so already plus much more
   */
  void init (PetscMatrix<T>* matrix);
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
   * Returns the raw PETSc preconditioner context pointer.  This allows
   * you to specify the PCShellSetApply() and PCShellSetSetUp() functions
   * if you desire.  Just don't do anything crazy like calling PCDestroy()!
   */
  PC pc() { this->init(); return _pc; }

  /**
   * Returns the raw PETSc ksp context pointer.  This is useful if
   * you are for example setting a custom convergence test with
   * KSPSetConvergenceTest().
   */
  KSP ksp() { this->init(); return _ksp; }
  
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
   * Tells PETSC to use the user-specified solver stored in
   * \p _solver_type
   */
  void set_petsc_solver_type ();

  /**
   * Internal function if shell matrix mode is used.
   */
  static PetscErrorCode _petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest);

  /**
   * Internal function if shell matrix mode is used.
   */
  static PetscErrorCode _petsc_shell_matrix_get_diagonal(Mat mat, Vec dest);

  // SLES removed from >= PETSc 2.2.0
#if PETSC_VERSION_LESS_THAN(2,2,0)
  
  /**
   * Linear solver context
   */
  SLES _sles;

#endif
  
  /**
   * Preconditioner context
   */
  PC _pc; 

  /**
   * Krylov subspace context
   */
  KSP _ksp;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
PetscLinearSolver<T>::PetscLinearSolver ()
{
  if (libMesh::n_processors() == 1)
    this->_preconditioner_type = ILU_PRECOND;
  else
    this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
}



template <typename T>
inline
PetscLinearSolver<T>::~PetscLinearSolver ()
{
  this->clear ();
}



#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
