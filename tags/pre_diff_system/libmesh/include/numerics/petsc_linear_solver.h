// $Id: petsc_linear_solver.h,v 1.3 2006-03-14 17:09:47 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

// C++ includes

// Local includes
#include "linear_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"


/**
 * Petsc include files.
 */
#ifdef HAVE_PETSC

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
}
#else
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
#endif



/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2005
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
  
private:

  /**
   * Tells PETSC to use the user-specified solver stored in
   * \p _solver_type
   */
  void set_petsc_solver_type ();

  /**
   * Tells PETSC to use the user-specified preconditioner stored in
   * \p _preconditioner_type
   */
  void set_petsc_preconditioner_type ();

  // SLES removed from >= PETSc 2.2.0
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
  
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



#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
