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



#ifndef LIBMESH_PETSC_LINEAR_SOLVER_H
#define LIBMESH_PETSC_LINEAR_SOLVER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_macro.h"

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

// Local includes
#include "libmesh/linear_solver.h"

// C++ includes
#include <cstddef>
#include <vector>

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

#if PETSC_RELEASE_LESS_THAN(3,0,1)
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


namespace libMesh
{

// forward declarations
template <typename T> class PetscMatrix;

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
  PetscLinearSolver (const libMesh::Parallel::Communicator &comm
		     LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

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
   * After calling this method, all successive solves will be
   * restricted to the given set of dofs, which must contain local
   * dofs on each processor only and not contain any duplicates.  This
   * mode can be disabled by calling this method with \p dofs being a
   * \p NULL pointer.
   */
  virtual void restrict_solve_to (const std::vector<unsigned int>* const dofs,
				  const SubsetSolveMode subset_solve_mode=SUBSET_ZERO);

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
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  std::pair<unsigned int, Real>
  adjoint_solve (SparseMatrix<T>  &matrix_in,
	 NumericVector<T> &solution_in,
	 NumericVector<T> &rhs_in,
	 const double tol,
		 const unsigned int m_its);


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
  static PetscErrorCode _petsc_shell_matrix_mult_add(Mat mat, Vec arg, Vec add, Vec dest);

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

  /**
   * PETSc index set containing the dofs on which to solve (\p NULL
   * means solve on all dofs).
   */
  IS _restrict_solve_to_is;

  /**
   * PETSc index set, complement to \p _restrict_solve_to_is.  This
   * will be created on demand by the method \p
   * _create_complement_is().
   */
  IS _restrict_solve_to_is_complement;

  /**
   * Internal method that returns the local size of \p
   * _restrict_solve_to_is.
   */
  size_t _restrict_solve_to_is_local_size(void)const;

  /**
   * Creates \p _restrict_solve_to_is_complement to contain all
   * indices that are local in \p vec_in, except those that are
   * contained in \p _restrict_solve_to_is.
   */
  void _create_complement_is (const NumericVector<T> &vec_in);

  /**
   * If restrict-solve-to-subset mode is active, this member decides
   * what happens with the dofs outside the subset.
   */
  SubsetSolveMode _subset_solve_mode;

};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
  PetscLinearSolver<T>::PetscLinearSolver (const libMesh::Parallel::Communicator &comm):
    LinearSolver<T>(comm),
    _restrict_solve_to_is(NULL),
    _restrict_solve_to_is_complement(NULL),
    _subset_solve_mode(SUBSET_ZERO)
{
  if (this->n_processors() == 1)
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



template <typename T>
inline size_t
PetscLinearSolver<T>::
_restrict_solve_to_is_local_size(void)const
{
  libmesh_assert(_restrict_solve_to_is);

  PetscInt s;
  int ierr = ISGetLocalSize(_restrict_solve_to_is,&s);
  LIBMESH_CHKERRABORT(ierr);

  return static_cast<size_t>(s);
}



template <typename T>
void
PetscLinearSolver<T>::_create_complement_is (const NumericVector<T> &
#if PETSC_VERSION_LESS_THAN(3,0,0)
                                             // unnamed to avoid compiler "unused parameter" warning
#else
                                             vec_in
#endif
  )
{
  libmesh_assert(_restrict_solve_to_is);
#if PETSC_VERSION_LESS_THAN(3,0,0)
  // No ISComplement in PETSc 2.3.3
  libmesh_not_implemented();
#else
  if(_restrict_solve_to_is_complement==NULL)
    {
      int ierr = ISComplement(_restrict_solve_to_is,
			      vec_in.first_local_index(),
			      vec_in.last_local_index(),
			      &_restrict_solve_to_is_complement);
      LIBMESH_CHKERRABORT(ierr);
    }
#endif
}



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_LINEAR_SOLVER_H
