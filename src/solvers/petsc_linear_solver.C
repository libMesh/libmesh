// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC


// C++ includes
#include <string.h>

// Local Includes
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/petsc_auto_fieldsplit.h"
#include "libmesh/solver_configuration.h"

namespace libMesh
{

extern "C"
{
#if PETSC_RELEASE_LESS_THAN(3,0,1)
  PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx)
  {
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number> *>(ctx);

    if (!preconditioner->initialized())
      libmesh_error_msg("Preconditioner not initialized!  Make sure you call init() before solve!");

    preconditioner->setup();

    return 0;
  }


  PetscErrorCode __libmesh_petsc_preconditioner_apply(void * ctx, Vec x, Vec y)
  {
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number> *>(ctx);

    PetscVector<Number> x_vec(x, preconditioner->comm());
    PetscVector<Number> y_vec(y, preconditioner->comm());

    preconditioner->apply(x_vec,y_vec);

    return 0;
  }
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC pc)
  {
    void * ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number> *>(ctx);

    if (!preconditioner->initialized())
      libmesh_error_msg("Preconditioner not initialized!  Make sure you call init() before solve!");

    preconditioner->setup();

    return 0;
  }

  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y)
  {
    void * ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number> *>(ctx);

    PetscVector<Number> x_vec(x, preconditioner->comm());
    PetscVector<Number> y_vec(y, preconditioner->comm());

    preconditioner->apply(x_vec,y_vec);

    return 0;
  }
#endif
} // end extern "C"

/*----------------------- functions ----------------------------------*/
template <typename T>
void PetscLinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      /* If we were restricted to some subset, this restriction must
         be removed and the subset index set destroyed.  */
      if (_restrict_solve_to_is != libmesh_nullptr)
        {
          PetscErrorCode ierr = LibMeshISDestroy(&_restrict_solve_to_is);
          LIBMESH_CHKERR(ierr);
          _restrict_solve_to_is = libmesh_nullptr;
        }

      if(_restrict_solve_to_is_complement != libmesh_nullptr)
        {
          PetscErrorCode ierr = LibMeshISDestroy(&_restrict_solve_to_is_complement);
          LIBMESH_CHKERR(ierr);
          _restrict_solve_to_is_complement = libmesh_nullptr;
        }

      this->_is_initialized = false;

      PetscErrorCode ierr=0;

      ierr = LibMeshKSPDestroy(&_ksp);
      LIBMESH_CHKERR(ierr);

      // Mimic PETSc default solver and preconditioner
      this->_solver_type           = GMRES;

      if(!this->_preconditioner)
        {
          if (this->n_processors() == 1)
            this->_preconditioner_type = ILU_PRECOND;
          else
            this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
        }
    }
}



template <typename T>
void PetscLinearSolver<T>::init (const char * name)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

      // Create the linear solver context
      ierr = KSPCreate (this->comm().get(), &_ksp);
      LIBMESH_CHKERR(ierr);

      if (name)
        {
          ierr = KSPSetOptionsPrefix(_ksp, name);
          LIBMESH_CHKERR(ierr);
        }

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
      LIBMESH_CHKERR(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // If the SolverConfiguration object is provided, use it to set
      // options during solver initialization.
      if(this->_solver_configuration)
        {
          this->_solver_configuration->set_options_during_init();
        }

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      LIBMESH_CHKERR(ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //LIBMESH_CHKERR(ierr);

      // Have the Krylov subspace method use our good initial guess
      // rather than 0, unless the user requested a KSPType of
      // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_RELEASE_LESS_THAN(3,4,0)
      // Pre-3.0 and petsc-dev (as of October 2012) use non-const versions
      KSPType ksp_type;
#else
      const KSPType ksp_type;
#endif

      ierr = KSPGetType (_ksp, &ksp_type);
      LIBMESH_CHKERR(ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
          LIBMESH_CHKERR(ierr);
        }

      // Notify PETSc of location to store residual history.
      // This needs to be called before any solves, since
      // it sets the residual history length to zero.  The default
      // behavior is for PETSc to allocate (internally) an array
      // of size 1000 to hold the residual norm history.
      ierr = KSPSetResidualHistory(_ksp,
                                   PETSC_NULL,   // pointer to the array which holds the history
                                   PETSC_DECIDE, // size of the array holding the history
                                   PETSC_TRUE);  // Whether or not to reset the history for each solve.
      LIBMESH_CHKERR(ierr);

      PetscPreconditioner<T>::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

      //If there is a preconditioner object we need to set the internal setup and apply routines
      if(this->_preconditioner)
        {
          this->_preconditioner->init();
          PCShellSetContext(_pc,(void *)this->_preconditioner);
          PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
          PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
        }
    }
}


template <typename T>
void PetscLinearSolver<T>::init (PetscMatrix<T> * matrix,
                                 const char * name)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

      // Create the linear solver context
      ierr = KSPCreate (this->comm().get(), &_ksp);
      LIBMESH_CHKERR(ierr);

      if (name)
        {
          ierr = KSPSetOptionsPrefix(_ksp, name);
          LIBMESH_CHKERR(ierr);
        }

      //ierr = PCCreate (this->comm().get(), &_pc);
      //     LIBMESH_CHKERR(ierr);

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
      LIBMESH_CHKERR(ierr);

      // Set operators. The input matrix works as the preconditioning matrix
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat());
#endif
      LIBMESH_CHKERR(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // If the SolverConfiguration object is provided, use it to set
      // options during solver initialization.
      if(this->_solver_configuration)
        {
          this->_solver_configuration->set_options_during_init();
        }

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      LIBMESH_CHKERR(ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //LIBMESH_CHKERR(ierr);

      // Have the Krylov subspace method use our good initial guess
      // rather than 0, unless the user requested a KSPType of
      // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_RELEASE_LESS_THAN(3,4,0)
      KSPType ksp_type;
#else
      const KSPType ksp_type;
#endif

      ierr = KSPGetType (_ksp, &ksp_type);
      LIBMESH_CHKERR(ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
          LIBMESH_CHKERR(ierr);
        }

      // Notify PETSc of location to store residual history.
      // This needs to be called before any solves, since
      // it sets the residual history length to zero.  The default
      // behavior is for PETSc to allocate (internally) an array
      // of size 1000 to hold the residual norm history.
      ierr = KSPSetResidualHistory(_ksp,
                                   PETSC_NULL,   // pointer to the array which holds the history
                                   PETSC_DECIDE, // size of the array holding the history
                                   PETSC_TRUE);  // Whether or not to reset the history for each solve.
      LIBMESH_CHKERR(ierr);

      PetscPreconditioner<T>::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
      if(this->_preconditioner)
        {
          this->_preconditioner->set_matrix(*matrix);
          this->_preconditioner->init();
          PCShellSetContext(_pc,(void *)this->_preconditioner);
          PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
          PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
        }
    }
}



template <typename T>
void
PetscLinearSolver<T>::init_names (const System & sys)
{
  petsc_auto_fieldsplit(this->pc(), sys);
}



template <typename T>
void
PetscLinearSolver<T>::restrict_solve_to (const std::vector<unsigned int> * const dofs,
                                         const SubsetSolveMode subset_solve_mode)
{
  PetscErrorCode ierr=0;

  /* The preconditioner (in particular if a default preconditioner)
     will have to be reset.  We call this->clear() to do that.  This
     call will also remove and free any previous subset that may have
     been set before.  */
  this->clear();

  _subset_solve_mode = subset_solve_mode;

  if (dofs != libmesh_nullptr)
    {
      PetscInt * petsc_dofs = libmesh_nullptr;
      ierr = PetscMalloc(dofs->size()*sizeof(PetscInt), &petsc_dofs);
      LIBMESH_CHKERR(ierr);

      for (std::size_t i=0; i<dofs->size(); i++)
        petsc_dofs[i] = (*dofs)[i];

      ierr = ISCreateLibMesh(this->comm().get(),dofs->size(),petsc_dofs,PETSC_OWN_POINTER,&_restrict_solve_to_is);
      LIBMESH_CHKERR(ierr);
    }
}



template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (SparseMatrix<T> &  matrix_in,
                             SparseMatrix<T> &  precond_in,
                             NumericVector<T> & solution_in,
                             NumericVector<T> & rhs_in,
                             const double tol,
                             const unsigned int m_its)
{
  LOG_SCOPE("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T> * matrix   = cast_ptr<PetscMatrix<T> *>(&matrix_in);
  PetscMatrix<T> * precond  = cast_ptr<PetscMatrix<T> *>(&precond_in);
  PetscVector<T> * solution = cast_ptr<PetscVector<T> *>(&solution_in);
  PetscVector<T> * rhs      = cast_ptr<PetscVector<T> *>(&rhs_in);

  this->init (matrix);

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close ();
  precond->close ();
  solution->close ();
  rhs->close ();

  //   // If matrix != precond, then this means we have specified a
  //   // special preconditioner, so reset preconditioner type to PCMAT.
  //   if (matrix != precond)
  //     {
  //       this->_preconditioner_type = USER_PRECOND;
  //       this->set_petsc_preconditioner_type ();
  //     }

  Mat submat = libmesh_nullptr;
  Mat subprecond = libmesh_nullptr;
  Vec subrhs = libmesh_nullptr;
  Vec subsolution = libmesh_nullptr;
  VecScatter scatter = libmesh_nullptr;
  UniquePtr<PetscMatrix<Number> > subprecond_matrix;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      PetscInt is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERR(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterCreate(rhs->vec(), _restrict_solve_to_is, subrhs, libmesh_nullptr, &scatter);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERR(ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERR(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
         system, we will now change the right hand side to compensate
         for this.  Note that this is not necessary if \p SUBSET_ZERO
         has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
        {
          _create_complement_is(rhs_in);
          PetscInt is_complement_local_size =
            cast_int<PetscInt>(rhs_in.local_size()-is_local_size);

          Vec subvec1 = libmesh_nullptr;
          Mat submat1 = libmesh_nullptr;
          VecScatter scatter1 = libmesh_nullptr;

          ierr = VecCreate(this->comm().get(),&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1, libmesh_nullptr, &scatter1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);
          ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);

          ierr = VecScale(subvec1,-1.0);
          LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
          ierr = MatGetSubMatrix(matrix->mat(),
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#else
          ierr = MatGetSubMatrix(matrix->mat(),
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#endif

          ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
          LIBMESH_CHKERR(ierr);

          ierr = LibMeshVecScatterDestroy(&scatter1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshVecDestroy(&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshMatDestroy(&submat1);
          LIBMESH_CHKERR(ierr);
        }
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, submat, subprecond,
                             this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, submat, subprecond);

      PetscBool ksp_reuse_preconditioner = this->same_preconditioner ? PETSC_TRUE : PETSC_FALSE;
      ierr = KSPSetReusePreconditioner(_ksp, ksp_reuse_preconditioner);
#endif
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          subprecond_matrix.reset(new PetscMatrix<Number>(subprecond, this->comm()));
          this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
        }
    }
  else
    {
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
                             this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat());

      PetscBool ksp_reuse_preconditioner = this->same_preconditioner ? PETSC_TRUE : PETSC_FALSE;
      ierr = KSPSetReusePreconditioner(_ksp, ksp_reuse_preconditioner);
#endif
      LIBMESH_CHKERR(ierr);

      if(this->_preconditioner)
        {
          this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
        }
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
                           PETSC_DEFAULT, max_its);
  LIBMESH_CHKERR(ierr);

  // Allow command line options to override anything set programmatically.
  ierr = KSPSetFromOptions(_ksp);
  LIBMESH_CHKERR(ierr);

  // If the SolverConfiguration object is provided, use it to override
  // solver options.
  if(this->_solver_configuration)
    {
      this->_solver_configuration->configure_solver();
    }

  // Solve the linear system
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERR(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERR(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERR(ierr);

  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      switch(_subset_solve_mode)
        {
        case SUBSET_ZERO:
          ierr = VecZeroEntries(solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_COPY_RHS:
          ierr = VecCopy(rhs->vec(),solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_DONT_TOUCH:
          /* Nothing to do here.  */
          break;

        default:
          libmesh_error_msg("Invalid subset solve mode = " << _subset_solve_mode);
        }
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          // Before subprecond_matrix gets cleaned up, we should give
          // the _preconditioner a different matrix.
          this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
        }

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERR(ierr);
    }

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}

template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::adjoint_solve (SparseMatrix<T> &  matrix_in,
                                     NumericVector<T> & solution_in,
                                     NumericVector<T> & rhs_in,
                                     const double tol,
                                     const unsigned int m_its)
{
  LOG_SCOPE("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T> * matrix   = cast_ptr<PetscMatrix<T> *>(&matrix_in);
  // Note that the matrix and precond matrix are the same
  PetscMatrix<T> * precond  = cast_ptr<PetscMatrix<T> *>(&matrix_in);
  PetscVector<T> * solution = cast_ptr<PetscVector<T> *>(&solution_in);
  PetscVector<T> * rhs      = cast_ptr<PetscVector<T> *>(&rhs_in);

  this->init (matrix);

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close ();
  precond->close ();
  solution->close ();
  rhs->close ();

  Mat submat = libmesh_nullptr;
  Mat subprecond = libmesh_nullptr;
  Vec subrhs = libmesh_nullptr;
  Vec subsolution = libmesh_nullptr;
  VecScatter scatter = libmesh_nullptr;
  UniquePtr<PetscMatrix<Number> > subprecond_matrix;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      PetscInt is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERR(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs, libmesh_nullptr, &scatter);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERR(ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERR(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
         system, we will now change the right hand side to compensate
         for this.  Note that this is not necessary if \p SUBSET_ZERO
         has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
        {
          _create_complement_is(rhs_in);
          PetscInt is_complement_local_size =
            cast_int<PetscInt>(rhs_in.local_size()-is_local_size);

          Vec subvec1 = libmesh_nullptr;
          Mat submat1 = libmesh_nullptr;
          VecScatter scatter1 = libmesh_nullptr;

          ierr = VecCreate(this->comm().get(),&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1, libmesh_nullptr, &scatter1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);
          ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);

          ierr = VecScale(subvec1,-1.0);
          LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
          ierr = MatGetSubMatrix(matrix->mat(),
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#else
          ierr = MatGetSubMatrix(matrix->mat(),
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#endif

          ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
          LIBMESH_CHKERR(ierr);

          ierr = LibMeshVecScatterDestroy(&scatter1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshVecDestroy(&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshMatDestroy(&submat1);
          LIBMESH_CHKERR(ierr);
        }
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, submat, subprecond,
                             this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, submat, subprecond);

      PetscBool ksp_reuse_preconditioner = this->same_preconditioner ? PETSC_TRUE : PETSC_FALSE;
      ierr = KSPSetReusePreconditioner(_ksp, ksp_reuse_preconditioner);
#endif
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          subprecond_matrix.reset(new PetscMatrix<Number>(subprecond, this->comm()));
          this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
        }
    }
  else
    {
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
                             this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat());

      PetscBool ksp_reuse_preconditioner = this->same_preconditioner ? PETSC_TRUE : PETSC_FALSE;
      ierr = KSPSetReusePreconditioner(_ksp, ksp_reuse_preconditioner);
#endif
      LIBMESH_CHKERR(ierr);

      if(this->_preconditioner)
        {
          this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
        }
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
                           PETSC_DEFAULT, max_its);
  LIBMESH_CHKERR(ierr);

  // Allow command line options to override anything set programmatically.
  ierr = KSPSetFromOptions(_ksp);
  LIBMESH_CHKERR(ierr);

  // Solve the linear system
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      ierr = KSPSolveTranspose (_ksp, subrhs, subsolution);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      ierr = KSPSolveTranspose (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERR(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERR(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERR(ierr);

  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      switch(_subset_solve_mode)
        {
        case SUBSET_ZERO:
          ierr = VecZeroEntries(solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_COPY_RHS:
          ierr = VecCopy(rhs->vec(),solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_DONT_TOUCH:
          /* Nothing to do here.  */
          break;

        default:
          libmesh_error_msg("Invalid subset solve mode = " << _subset_solve_mode);
        }
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          // Before subprecond_matrix gets cleaned up, we should give
          // the _preconditioner a different matrix.
          this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
        }

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERR(ierr);
    }

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}


template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (const ShellMatrix<T> & shell_matrix,
                             NumericVector<T> & solution_in,
                             NumericVector<T> & rhs_in,
                             const double tol,
                             const unsigned int m_its)
{

#if PETSC_VERSION_LESS_THAN(3,1,0)
  if (_restrict_solve_to_is != libmesh_nullptr)
    libmesh_error_msg("The current implementation of subset solves with " \
                      << "shell matrices requires PETSc version 3.1 or above.  " \
                      << "Older PETSc version do not support automatic " \
                      << "submatrix generation of shell matrices.");
#endif

  LOG_SCOPE("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscVector<T> * solution = cast_ptr<PetscVector<T> *>(&solution_in);
  PetscVector<T> * rhs      = cast_ptr<PetscVector<T> *>(&rhs_in);

  this->init ();

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  Mat submat = libmesh_nullptr;
  Vec subrhs = libmesh_nullptr;
  Vec subsolution = libmesh_nullptr;
  VecScatter scatter = libmesh_nullptr;

  // Close the matrices and vectors in case this wasn't already done.
  solution->close ();
  rhs->close ();

  // Prepare the matrix.
  Mat mat;
  ierr = MatCreateShell(this->comm().get(),
                        rhs_in.local_size(),
                        solution_in.local_size(),
                        rhs_in.size(),
                        solution_in.size(),
                        const_cast<void *>(static_cast<const void *>(&shell_matrix)),
                        &mat);
  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void *.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T> *.  */

  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      PetscInt is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERR(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs, libmesh_nullptr, &scatter);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
         system, we will now change the right hand side to compensate
         for this.  Note that this is not necessary if \p SUBSET_ZERO
         has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
        {
          _create_complement_is(rhs_in);
          PetscInt is_complement_local_size =
            cast_int<PetscInt>(rhs_in.local_size()-is_local_size);

          Vec subvec1 = libmesh_nullptr;
          Mat submat1 = libmesh_nullptr;
          VecScatter scatter1 = libmesh_nullptr;

          ierr = VecCreate(this->comm().get(),&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1, libmesh_nullptr, &scatter1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);
          ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);

          ierr = VecScale(subvec1,-1.0);
          LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
          /* This point can't be reached, see above.  */
          libmesh_assert(false);
#else
          ierr = MatGetSubMatrix(mat,
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#endif

          // The following lines would be correct, but don't work
          // correctly in PETSc up to 3.1.0-p5.  See discussion in
          // petsc-users of Nov 9, 2010.
          //
          // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
          // LIBMESH_CHKERR(ierr);
          //
          // We workaround by using a temporary vector.  Note that the
          // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
          // so this is no effective performance loss.
          Vec subvec2 = libmesh_nullptr;
          ierr = VecCreate(this->comm().get(),&subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = MatMult(submat1,subvec1,subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = VecAXPY(subrhs,1.0,subvec2);

          ierr = LibMeshVecScatterDestroy(&scatter1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshVecDestroy(&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshMatDestroy(&submat1);
          LIBMESH_CHKERR(ierr);
        }
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, submat, submat,
                             DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, submat, submat);
#endif
      LIBMESH_CHKERR(ierr);
    }
  else
    {
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, mat, mat,
                             DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, mat, mat);
#endif
      LIBMESH_CHKERR(ierr);
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
                           PETSC_DEFAULT, max_its);
  LIBMESH_CHKERR(ierr);

  // Allow command line options to override anything set programmatically.
  ierr = KSPSetFromOptions(_ksp);
  LIBMESH_CHKERR(ierr);

  // Solve the linear system
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERR(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERR(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERR(ierr);

  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      switch(_subset_solve_mode)
        {
        case SUBSET_ZERO:
          ierr = VecZeroEntries(solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_COPY_RHS:
          ierr = VecCopy(rhs->vec(),solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_DONT_TOUCH:
          /* Nothing to do here.  */
          break;

        default:
          libmesh_error_msg("Invalid subset solve mode = " << _subset_solve_mode);
        }
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERR(ierr);

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERR(ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  LIBMESH_CHKERR(ierr);

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}



template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (const ShellMatrix<T> & shell_matrix,
                             const SparseMatrix<T> & precond_matrix,
                             NumericVector<T> & solution_in,
                             NumericVector<T> & rhs_in,
                             const double tol,
                             const unsigned int m_its)
{

#if PETSC_VERSION_LESS_THAN(3,1,0)
  if (_restrict_solve_to_is != libmesh_nullptr)
    libmesh_error_msg("The current implementation of subset solves with " \
                      << "shell matrices requires PETSc version 3.1 or above.  " \
                      << "Older PETSc version do not support automatic " \
                      << "submatrix generation of shell matrices.");
#endif

  LOG_SCOPE("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  const PetscMatrix<T> * precond  = cast_ptr<const PetscMatrix<T> *>(&precond_matrix);
  PetscVector<T> * solution = cast_ptr<PetscVector<T> *>(&solution_in);
  PetscVector<T> * rhs      = cast_ptr<PetscVector<T> *>(&rhs_in);

  this->init ();

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  Mat submat = libmesh_nullptr;
  Mat subprecond = libmesh_nullptr;
  Vec subrhs = libmesh_nullptr;
  Vec subsolution = libmesh_nullptr;
  VecScatter scatter = libmesh_nullptr;
  UniquePtr<PetscMatrix<Number> > subprecond_matrix;

  // Close the matrices and vectors in case this wasn't already done.
  solution->close ();
  rhs->close ();

  // Prepare the matrix.
  Mat mat;
  ierr = MatCreateShell(this->comm().get(),
                        rhs_in.local_size(),
                        solution_in.local_size(),
                        rhs_in.size(),
                        solution_in.size(),
                        const_cast<void *>(static_cast<const void *>(&shell_matrix)),
                        &mat);
  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void *.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T> *.  */

  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      PetscInt is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERR(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs, libmesh_nullptr, &scatter);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERR(ierr);
      ierr = MatGetSubMatrix(const_cast<PetscMatrix<T> *>(precond)->mat(),
                             _restrict_solve_to_is,_restrict_solve_to_is,
                             MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERR(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
         system, we will now change the right hand side to compensate
         for this.  Note that this is not necessary if \p SUBSET_ZERO
         has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
        {
          _create_complement_is(rhs_in);
          PetscInt is_complement_local_size = rhs_in.local_size()-is_local_size;

          Vec subvec1 = libmesh_nullptr;
          Mat submat1 = libmesh_nullptr;
          VecScatter scatter1 = libmesh_nullptr;

          ierr = VecCreate(this->comm().get(),&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1, libmesh_nullptr, &scatter1);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);
          ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);

          ierr = VecScale(subvec1,-1.0);
          LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
          /* This point can't be reached, see above.  */
          libmesh_assert(false);
#else
          ierr = MatGetSubMatrix(mat,
                                 _restrict_solve_to_is,_restrict_solve_to_is_complement,
                                 MAT_INITIAL_MATRIX,&submat1);
          LIBMESH_CHKERR(ierr);
#endif

          // The following lines would be correct, but don't work
          // correctly in PETSc up to 3.1.0-p5.  See discussion in
          // petsc-users of Nov 9, 2010.
          //
          // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
          // LIBMESH_CHKERR(ierr);
          //
          // We workaround by using a temporary vector.  Note that the
          // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
          // so this is no effective performance loss.
          Vec subvec2 = libmesh_nullptr;
          ierr = VecCreate(this->comm().get(),&subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
          LIBMESH_CHKERR(ierr);
          ierr = VecSetFromOptions(subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = MatMult(submat1,subvec1,subvec2);
          LIBMESH_CHKERR(ierr);
          ierr = VecAXPY(subrhs,1.0,subvec2);
          LIBMESH_CHKERR(ierr);

          ierr = LibMeshVecScatterDestroy(&scatter1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshVecDestroy(&subvec1);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshMatDestroy(&submat1);
          LIBMESH_CHKERR(ierr);
        }

#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, submat, subprecond,
                             DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, submat, subprecond);
#endif
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          subprecond_matrix.reset(new PetscMatrix<Number>(subprecond, this->comm()));
          this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
        }
    }
  else
    {
#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix<T> *>(precond)->mat(),
                             DIFFERENT_NONZERO_PATTERN);
#else
      ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix<T> *>(precond)->mat());
#endif
      LIBMESH_CHKERR(ierr);

      if(this->_preconditioner)
        {
          this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number> &>(precond_matrix));
          this->_preconditioner->init();
        }
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
                           PETSC_DEFAULT, max_its);
  LIBMESH_CHKERR(ierr);

  // Allow command line options to override anything set programmatically.
  ierr = KSPSetFromOptions(_ksp);
  LIBMESH_CHKERR(ierr);

  // Solve the linear system
  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERR(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERR(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERR(ierr);

  if (_restrict_solve_to_is != libmesh_nullptr)
    {
      switch(_subset_solve_mode)
        {
        case SUBSET_ZERO:
          ierr = VecZeroEntries(solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_COPY_RHS:
          ierr = VecCopy(rhs->vec(),solution->vec());
          LIBMESH_CHKERR(ierr);
          break;

        case SUBSET_DONT_TOUCH:
          /* Nothing to do here.  */
          break;

        default:
          libmesh_error_msg("Invalid subset solve mode = " << _subset_solve_mode);
        }
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERR(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERR(ierr);

      if (this->_preconditioner)
        {
          // Before subprecond_matrix gets cleaned up, we should give
          // the _preconditioner a different matrix.
          this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number> &>(precond_matrix));
          this->_preconditioner->init();
        }

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERR(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERR(ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  LIBMESH_CHKERR(ierr);

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}



template <typename T>
void PetscLinearSolver<T>::get_residual_history(std::vector<double> & hist)
{
  PetscErrorCode ierr = 0;
  PetscInt its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal * p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  LIBMESH_CHKERR(ierr);

  // Check for early return
  if (its == 0) return;

  // Create space to store the result
  hist.resize(its);

  // Copy history into the vector provided by the user.
  for (PetscInt i=0; i<its; ++i)
    {
      hist[i] = *p;
      p++;
    }
}




template <typename T>
Real PetscLinearSolver<T>::get_initial_residual()
{
  PetscErrorCode ierr = 0;
  PetscInt its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal * p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  LIBMESH_CHKERR(ierr);

  // Check no residual history
  if (its == 0)
    {
      libMesh::err << "No iterations have been performed, returning 0." << std::endl;
      return 0.;
    }

  // Otherwise, return the value pointed to by p.
  return *p;
}




template <typename T>
void PetscLinearSolver<T>::set_petsc_solver_type()
{
  PetscErrorCode ierr = 0;

  switch (this->_solver_type)
    {

    case CG:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPCG));
      LIBMESH_CHKERR(ierr);
      return;

    case CR:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPCR));
      LIBMESH_CHKERR(ierr);
      return;

    case CGS:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPCGS));
      LIBMESH_CHKERR(ierr);
      return;

    case BICG:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPBICG));
      LIBMESH_CHKERR(ierr);
      return;

    case TCQMR:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPTCQMR));
      LIBMESH_CHKERR(ierr);
      return;

    case TFQMR:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPTFQMR));
      LIBMESH_CHKERR(ierr);
      return;

    case LSQR:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPLSQR));
      LIBMESH_CHKERR(ierr);
      return;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPBCGS));
      LIBMESH_CHKERR(ierr);
      return;

    case MINRES:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPMINRES));
      LIBMESH_CHKERR(ierr);
      return;

    case GMRES:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPGMRES));
      LIBMESH_CHKERR(ierr);
      return;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPRICHARDSON));
      LIBMESH_CHKERR(ierr);
      return;

    case CHEBYSHEV:
#if defined(LIBMESH_HAVE_PETSC) && PETSC_VERSION_LESS_THAN(3,3,0)
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPCHEBYCHEV));
      LIBMESH_CHKERR(ierr);
      return;
#else
      ierr = KSPSetType (_ksp, const_cast<KSPType>(KSPCHEBYSHEV));
      LIBMESH_CHKERR(ierr);
      return;
#endif


    default:
      libMesh::err << "ERROR:  Unsupported PETSC Solver: "
                   << Utility::enum_to_string(this->_solver_type) << std::endl
                   << "Continuing with PETSC defaults" << std::endl;
    }
}



template <typename T>
LinearConvergenceReason PetscLinearSolver<T>::get_converged_reason() const
{
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);

  switch(reason)
    {
#if !PETSC_VERSION_LESS_THAN(3,2,0)
    case KSP_CONVERGED_RTOL_NORMAL     : return CONVERGED_RTOL_NORMAL;
    case KSP_CONVERGED_ATOL_NORMAL     : return CONVERGED_ATOL_NORMAL;
#endif
    case KSP_CONVERGED_RTOL            : return CONVERGED_RTOL;
    case KSP_CONVERGED_ATOL            : return CONVERGED_ATOL;
    case KSP_CONVERGED_ITS             : return CONVERGED_ITS;
    case KSP_CONVERGED_CG_NEG_CURVE    : return CONVERGED_CG_NEG_CURVE;
    case KSP_CONVERGED_CG_CONSTRAINED  : return CONVERGED_CG_CONSTRAINED;
    case KSP_CONVERGED_STEP_LENGTH     : return CONVERGED_STEP_LENGTH;
    case KSP_CONVERGED_HAPPY_BREAKDOWN : return CONVERGED_HAPPY_BREAKDOWN;
    case KSP_DIVERGED_NULL             : return DIVERGED_NULL;
    case KSP_DIVERGED_ITS              : return DIVERGED_ITS;
    case KSP_DIVERGED_DTOL             : return DIVERGED_DTOL;
    case KSP_DIVERGED_BREAKDOWN        : return DIVERGED_BREAKDOWN;
    case KSP_DIVERGED_BREAKDOWN_BICG   : return DIVERGED_BREAKDOWN_BICG;
    case KSP_DIVERGED_NONSYMMETRIC     : return DIVERGED_NONSYMMETRIC;
    case KSP_DIVERGED_INDEFINITE_PC    : return DIVERGED_INDEFINITE_PC;
#if PETSC_VERSION_LESS_THAN(3,4,0)
    case KSP_DIVERGED_NAN              : return DIVERGED_NAN;
#else
    case KSP_DIVERGED_NANORINF         : return DIVERGED_NAN;
#endif
    case KSP_DIVERGED_INDEFINITE_MAT   : return DIVERGED_INDEFINITE_MAT;
    case KSP_CONVERGED_ITERATING       : return CONVERGED_ITERATING;
#if !PETSC_VERSION_LESS_THAN(3,7,0)
    case KSP_DIVERGED_PCSETUP_FAILED   : return DIVERGED_PCSETUP_FAILED;
#endif
    default :
      libMesh::err << "Unknown convergence flag!" << std::endl;
      return UNKNOWN_FLAG;
    }
}


template <typename T>
PetscErrorCode PetscLinearSolver<T>::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
{
  /* Get the matrix context.  */
  PetscErrorCode ierr=0;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);
  CHKERRABORT(shell_matrix.comm().get(), ierr);

  /* Make \p NumericVector instances around the vectors.  */
  PetscVector<T> arg_global(arg, shell_matrix.comm());
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  /* Call the user function.  */
  shell_matrix.vector_mult(dest_global,arg_global);

  return ierr;
}



template <typename T>
PetscErrorCode PetscLinearSolver<T>::_petsc_shell_matrix_mult_add(Mat mat, Vec arg, Vec add, Vec dest)
{
  /* Get the matrix context.  */
  PetscErrorCode ierr=0;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);
  CHKERRABORT(shell_matrix.comm().get(), ierr);

  /* Make \p NumericVector instances around the vectors.  */
  PetscVector<T> arg_global(arg, shell_matrix.comm());
  PetscVector<T> dest_global(dest, shell_matrix.comm());
  PetscVector<T> add_global(add, shell_matrix.comm());

  if(add!=arg)
    {
      arg_global = add_global;
    }

  /* Call the user function.  */
  shell_matrix.vector_mult_add(dest_global,arg_global);

  return ierr;
}



template <typename T>
PetscErrorCode PetscLinearSolver<T>::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest)
{
  /* Get the matrix context.  */
  PetscErrorCode ierr=0;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);
  CHKERRABORT(shell_matrix.comm().get(), ierr);

  /* Make \p NumericVector instances around the vector.  */
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  /* Call the user function.  */
  shell_matrix.get_diagonal(dest_global);

  return ierr;
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscLinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
