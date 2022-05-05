// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/wrapped_petsc.h"

// Petsc include files.
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscksp.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

// Local includes
#include "libmesh/linear_solver.h"

// C++ includes
#include <cstddef>
#include <vector>
#include <memory>

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed for preconditioning
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  /**
   * This function is called by PETSc to initialize the preconditioner.
   * ctx will hold the Preconditioner.
   */
  PetscErrorCode libmesh_petsc_preconditioner_setup (PC);

  /**
   * This function is called by PETSc to actually apply the preconditioner.
   * ctx will hold the Preconditioner.
   */
  PetscErrorCode libmesh_petsc_preconditioner_apply(PC, Vec x, Vec y);

#if LIBMESH_ENABLE_DEPRECATED
  /**
   * This function is called by PETSc to initialize the preconditioner.
   * ctx will hold the Preconditioner.
   */
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC);

  /**
   * This function is called by PETSc to actually apply the preconditioner.
   * ctx will hold the Preconditioner.
   */
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
 * \author Benjamin Kirk
 * \date 2002-2007
 */
template <typename T>
class PetscLinearSolver : public LinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Petsc data structures
   */
  PetscLinearSolver (const libMesh::Parallel::Communicator & comm_in);

  /**
   * Destructor.
   */
  virtual ~PetscLinearSolver () = default;

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () override;

  /**
   * Initialize data structures if not done so already.
   * Assigns a name, which is turned into an underscore-separated
   * prefix for the underlying KSP object.
   */
  virtual void init (const char * name = nullptr) override;

  /**
   * Initialize data structures if not done so already plus much more
   */
  void init (PetscMatrix<T> * matrix,
             const char * name = nullptr);

  /**
   * Apply names to the system to be solved.  This sets an option
   * prefix from the system name and sets field names from the
   * system's variable names.
   *
   * Since field names are applied to DoF numberings, this method must
   * be called again after any System reinit.
   */
  virtual void init_names (const System &) override;

  /**
   * After calling this method, all successive solves will be
   * restricted to the given set of dofs, which must contain local
   * dofs on each processor only and not contain any duplicates.  This
   * mode can be disabled by calling this method with \p dofs being a
   * \p nullptr.
   */
  virtual void restrict_solve_to (const std::vector<unsigned int> * const dofs,
                                  const SubsetSolveMode subset_solve_mode=SUBSET_ZERO) override;

  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real>
  solve (SparseMatrix<T> & matrix_in,
         NumericVector<T> & solution_in,
         NumericVector<T> & rhs_in,
         const double tol,
         const unsigned int m_its) override
  {
    return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
  }


  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real>
  adjoint_solve (SparseMatrix<T> & matrix_in,
                 NumericVector<T> & solution_in,
                 NumericVector<T> & rhs_in,
                 const double tol,
                 const unsigned int m_its) override;

  /**
   * This method allows you to call a linear solver while specifying
   * the matrix to use as the (left) preconditioning matrix.
   *
   * \note The linear solver will not compute a preconditioner in this
   * case, and will instead premultiply by the matrix you provide.  In
   * PETSc, this is accomplished by calling
   * \code
   * PCSetType(_pc, PCMAT);
   * \endcode
   * before invoking KSPSolve().
   *
   * \note This functionality is not implemented in the LinearSolver
   * class since there is not a built-in analog to this method for
   * LASPACK. You could probably implement it by hand if you wanted.
   */
  virtual std::pair<unsigned int, Real>
  solve (SparseMatrix<T> & matrix,
         SparseMatrix<T> & preconditioner,
         NumericVector<T> & solution,
         NumericVector<T> & rhs,
         const double tol,
         const unsigned int m_its) override;

  /**
   * This function solves a system whose matrix is a shell matrix.
   */
  virtual std::pair<unsigned int, Real>
  solve (const ShellMatrix<T> & shell_matrix,
         NumericVector<T> & solution_in,
         NumericVector<T> & rhs_in,
         const double tol,
         const unsigned int m_its) override;

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
         const unsigned int m_its) override;

  /**
   * \returns The raw PETSc preconditioner context pointer.
   *
   * This allows you to specify the PCShellSetApply() and
   * PCShellSetSetUp() functions if you desire.  Just don't do
   * anything crazy like calling libMeshPCDestroy() on the pointer!
   */
  PC pc() { this->init(); return _pc; }

  /**
   * \returns The raw PETSc ksp context pointer.
   *
   * This is useful if you are for example setting a custom
   * convergence test with KSPSetConvergenceTest().
   */
  KSP ksp();

  /**
   * Fills the input vector with the sequence of residual norms
   * from the latest iterative solve.
   */
  void get_residual_history(std::vector<double> & hist);

  /**
   * \returns Just the initial residual for the solve just
   * completed with this interface.
   *
   * Use this method instead of the one above if you just want the
   * starting residual and not the entire history.
   */
  Real get_initial_residual();

  /**
   * \returns The solver's convergence flag
   */
  virtual LinearConvergenceReason get_converged_reason() const override;

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

  /**
   * Preconditioner context
   */
  PC _pc;

  /**
   * Krylov subspace context
   */
  WrappedPetsc<KSP> _ksp;

  /**
   * PETSc index set containing the dofs on which to solve (\p nullptr
   * means solve on all dofs).
   */
  WrappedPetsc<IS> _restrict_solve_to_is;

  /**
   * PETSc index set, complement to \p _restrict_solve_to_is.  This
   * will be created on demand by the method \p
   * _create_complement_is().
   */
  WrappedPetsc<IS> _restrict_solve_to_is_complement;

  /**
   * \returns The local size of \p _restrict_solve_to_is.
   */
  PetscInt restrict_solve_to_is_local_size() const;

  /**
   * Creates \p _restrict_solve_to_is_complement to contain all
   * indices that are local in \p vec_in, except those that are
   * contained in \p _restrict_solve_to_is.
   */
  void create_complement_is (const NumericVector<T> & vec_in);

  /**
   * If restrict-solve-to-subset mode is active, this member decides
   * what happens with the dofs outside the subset.
   */
  SubsetSolveMode _subset_solve_mode;
};

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_LINEAR_SOLVER_H
