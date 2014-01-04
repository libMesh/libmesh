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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC


// C++ includes
#include <string.h>

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

extern "C"
{
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif


#if PETSC_RELEASE_LESS_THAN(3,0,1)
  PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx)
  {
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

    if(!preconditioner->initialized())
    {
      err<<"Preconditioner not initialized!  Make sure you call init() before solve!"<<std::endl;
      libmesh_error();
    }

    preconditioner->setup();

    return 0;
  }


  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)
  {
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

    PetscVector<Number> x_vec(x, preconditioner->comm());
    PetscVector<Number> y_vec(y, preconditioner->comm());

    preconditioner->apply(x_vec,y_vec);

    return 0;
  }
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC pc)
  {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

    if(!preconditioner->initialized())
    {
      err<<"Preconditioner not initialized!  Make sure you call init() before solve!"<<std::endl;
      libmesh_error();
    }

    preconditioner->setup();

    return 0;
  }

  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y)
  {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

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
      if(_restrict_solve_to_is!=NULL)
	{
	  PetscErrorCode ierr = LibMeshISDestroy(&_restrict_solve_to_is);
	  LIBMESH_CHKERRABORT(ierr);
	  _restrict_solve_to_is = NULL;
	}

      if(_restrict_solve_to_is_complement!=NULL)
	{
	  PetscErrorCode ierr = LibMeshISDestroy(&_restrict_solve_to_is_complement);
	  LIBMESH_CHKERRABORT(ierr);
	  _restrict_solve_to_is_complement = NULL;
	}

      this->_is_initialized = false;

      PetscErrorCode ierr=0;

#if PETSC_VERSION_LESS_THAN(2,2,0)

  // 2.1.x & earlier style
  ierr = SLESDestroy(_sles);
         LIBMESH_CHKERRABORT(ierr);

#else

  // 2.2.0 & newer style
  ierr = LibMeshKSPDestroy(&_ksp);
         LIBMESH_CHKERRABORT(ierr);

#endif


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
void PetscLinearSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

      // Create the linear solver context
      ierr = SLESCreate (this->comm().get(), &_sles);
             LIBMESH_CHKERRABORT(ierr);

      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp);
             LIBMESH_CHKERRABORT(ierr);
      ierr = SLESGetPC        (_sles, &_pc);
             LIBMESH_CHKERRABORT(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = SLESSetFromOptions (_sles);
             LIBMESH_CHKERRABORT(ierr);

// 2.2.0 & newer style
#else

      // Create the linear solver context
      ierr = KSPCreate (this->comm().get(), &_ksp);
             LIBMESH_CHKERRABORT(ierr);

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
             LIBMESH_CHKERRABORT(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      LIBMESH_CHKERRABORT(ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //LIBMESH_CHKERRABORT(ierr);


#endif

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
             LIBMESH_CHKERRABORT(ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
                 LIBMESH_CHKERRABORT(ierr);
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
      LIBMESH_CHKERRABORT(ierr);

      PetscPreconditioner<T>::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

      //If there is a preconditioner object we need to set the internal setup and apply routines
      if(this->_preconditioner)
      {
        this->_preconditioner->init();
        PCShellSetContext(_pc,(void*)this->_preconditioner);
        PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
        PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
      }
    }
}


template <typename T>
void PetscLinearSolver<T>::init ( PetscMatrix<T>* matrix )
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

      // Create the linear solver context
      ierr = SLESCreate (this->comm().get(), &_sles);
             LIBMESH_CHKERRABORT(ierr);

      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp);
             LIBMESH_CHKERRABORT(ierr);
      ierr = SLESGetPC        (_sles, &_pc);
             LIBMESH_CHKERRABORT(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = SLESSetFromOptions (_sles);
             LIBMESH_CHKERRABORT(ierr);

// 2.2.0 & newer style
#else

      // Create the linear solver context
      ierr = KSPCreate (this->comm().get(), &_ksp);
             LIBMESH_CHKERRABORT(ierr);


      //ierr = PCCreate (this->comm().get(), &_pc);
      //     LIBMESH_CHKERRABORT(ierr);

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
             LIBMESH_CHKERRABORT(ierr);

      // Set operators. The input matrix works as the preconditioning matrix
      ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),DIFFERENT_NONZERO_PATTERN);
             LIBMESH_CHKERRABORT(ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      LIBMESH_CHKERRABORT(ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //LIBMESH_CHKERRABORT(ierr);


#endif

      // Have the Krylov subspace method use our good initial guess
      // rather than 0, unless the user requested a KSPType of
      // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_RELEASE_LESS_THAN(3,4,0)
      KSPType ksp_type;
#else
      const KSPType ksp_type;
#endif

      ierr = KSPGetType (_ksp, &ksp_type);
             LIBMESH_CHKERRABORT(ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
                 LIBMESH_CHKERRABORT(ierr);
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
      LIBMESH_CHKERRABORT(ierr);

      PetscPreconditioner<T>::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
      if(this->_preconditioner)
      {
        this->_preconditioner->set_matrix(*matrix);
        this->_preconditioner->init();
        PCShellSetContext(_pc,(void*)this->_preconditioner);
        PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
        PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
      }
    }
}



template <typename T>
void
PetscLinearSolver<T>::restrict_solve_to (const std::vector<unsigned int>* const dofs,
					 const SubsetSolveMode subset_solve_mode)
{
  PetscErrorCode ierr=0;

  /* The preconditioner (in particular if a default preconditioner)
     will have to be reset.  We call this->clear() to do that.  This
     call will also remove and free any previous subset that may have
     been set before.  */
  this->clear();

  _subset_solve_mode = subset_solve_mode;

  if(dofs!=NULL)
    {
      PetscInt* petsc_dofs = NULL;
      ierr = PetscMalloc(dofs->size()*sizeof(PetscInt), &petsc_dofs);
      LIBMESH_CHKERRABORT(ierr);

      for(size_t i=0; i<dofs->size(); i++)
	{
	  petsc_dofs[i] = (*dofs)[i];
	}

      ierr = ISCreateLibMesh(this->comm().get(),dofs->size(),petsc_dofs,PETSC_OWN_POINTER,&_restrict_solve_to_is);
      LIBMESH_CHKERRABORT(ierr);
    }
}



template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (SparseMatrix<T>&  matrix_in,
			     SparseMatrix<T>&  precond_in,
			     NumericVector<T>& solution_in,
			     NumericVector<T>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{
  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T>* matrix   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_in);
  PetscMatrix<T>* precond  = libmesh_cast_ptr<PetscMatrix<T>*>(&precond_in);
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

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

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

  if(_restrict_solve_to_is!=NULL)
    {
      libmesh_not_implemented();
    }

  // Set operators. The input matrix works as the preconditioning matrix
  ierr = SLESSetOperators(_sles, matrix->mat(), precond->mat(),
			  DIFFERENT_NONZERO_PATTERN);
         LIBMESH_CHKERRABORT(ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
         LIBMESH_CHKERRABORT(ierr);


  // Solve the linear system
  ierr = SLESSolve (_sles, rhs->vec(), solution->vec(), &its);
         LIBMESH_CHKERRABORT(ierr);


  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         LIBMESH_CHKERRABORT(ierr);

// 2.2.0
#elif PETSC_VERSION_LESS_THAN(2,2,1)

  if(_restrict_solve_to_is!=NULL)
    {
      libmesh_not_implemented();
    }

  // Set operators. The input matrix works as the preconditioning matrix
  //ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
	 //			 SAME_NONZERO_PATTERN);
	 //       LIBMESH_CHKERRABORT(ierr);


  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  // Convergence is detected at iteration k if
  // ||r_k||_2 < max(rtol*||b||_2 , abstol)
  // where r_k is the residual vector and b is the right-hand side.  Note that
  // it is the *maximum* of the two values, the larger of which will almost
  // always be rtol*||b||_2.
  ierr = KSPSetTolerances (_ksp,
			   tol,           // rtol   = relative decrease in residual  (1.e-5)
			   PETSC_DEFAULT, // abstol = absolute convergence tolerance (1.e-50)
 			   PETSC_DEFAULT, // dtol   = divergence tolerance           (1.e+5)
			   max_its);
         LIBMESH_CHKERRABORT(ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec());
         LIBMESH_CHKERRABORT(ierr);

  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec());
         LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  ierr = KSPSolve (_ksp);
         LIBMESH_CHKERRABORT(ierr);

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         LIBMESH_CHKERRABORT(ierr);

// 2.2.1 & newer style
#else

  Mat submat = NULL;
  Mat subprecond = NULL;
  Vec subrhs = NULL;
  Vec subsolution = NULL;
  VecScatter scatter = NULL;
  PetscMatrix<Number>* subprecond_matrix = NULL;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERRABORT(ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERRABORT(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
	{
	  _create_complement_is(rhs_in);
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;
	  Mat submat1 = NULL;
	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScale(subvec1,-1.0);
	  LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#else
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#endif

	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  LIBMESH_CHKERRABORT(ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond,
						      this->comm());
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

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
         LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERRABORT(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERRABORT(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         LIBMESH_CHKERRABORT(ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  /* Before we delete subprecond_matrix, we should give the
	     _preconditioner a different matrix.  */
	  this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
	  delete subprecond_matrix;
	  subprecond_matrix = NULL;
	}

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERRABORT(ierr);
    }

#endif

  STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}

template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::adjoint_solve (SparseMatrix<T>&  matrix_in,
				     NumericVector<T>& solution_in,
				     NumericVector<T>& rhs_in,
				     const double tol,
				     const unsigned int m_its)
{
  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T>* matrix   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_in);
  // Note that the matrix and precond matrix are the same
  PetscMatrix<T>* precond  = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_in);
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

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

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

  if(_restrict_solve_to_is!=NULL)
    {
      libmesh_not_implemented();
    }

  // Based on http://wolfgang.math.tamu.edu/svn/public/deal.II/branches/MATH676/2008/deal.II/lac/source/petsc_solver.cc, http://tccc.iesl.forth.gr/AMS_EPEAEK/Elements/doc/in_html/petsc/SLES/index.html

  SLES sles;
  ierr = SLESCreate (this->comm().get(), &sles);
  LIBMESH_CHKERRABORT(ierr);

  ierr = SLESSetOperators (sles, matrix->mat(), precond->mat(), this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
  LIBMESH_CHKERRABORT(ierr);

  KSP ksp;
  ierr = SLESGetKSP (sles, &ksp);
  LIBMESH_CHKERRABORT(ierr);

  ierr = SLESSetUp (sles, rhs->vec(), solution->vec());
  LIBMESH_CHKERRABORT(ierr);

  // See http://tccc.iesl.forth.gr/AMS_EPEAEK/Elements/doc/in_html/petsc/KSP/KSPSolveTrans.html#KSPSolveTrans
  ierr = SLESSolveTrans (ksp, &its);
  LIBMESH_CHKERRABORT(ierr);

// 2.2.0
#elif PETSC_VERSION_LESS_THAN(2,2,1)

  if(_restrict_solve_to_is!=NULL)
    {
      libmesh_not_implemented();
    }

  // Set operators. The input matrix works as the preconditioning matrix
  // This was commented earlier but it looks like KSPSetOperators is supported
  // after PETSc 2.2.0
  ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			 this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
  LIBMESH_CHKERRABORT(ierr);


  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  // Convergence is detected at iteration k if
  // ||r_k||_2 < max(rtol*||b||_2 , abstol)
  // where r_k is the residual vector and b is the right-hand side.  Note that
  // it is the *maximum* of the two values, the larger of which will almost
  // always be rtol*||b||_2.
  ierr = KSPSetTolerances (_ksp,
			   tol,           // rtol   = relative decrease in residual  (1.e-5)
			   PETSC_DEFAULT, // abstol = absolute convergence tolerance (1.e-50)
 			   PETSC_DEFAULT, // dtol   = divergence tolerance           (1.e+5)
			   max_its);
         LIBMESH_CHKERRABORT(ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec());
         LIBMESH_CHKERRABORT(ierr);

  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec());
         LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  ierr = KSPSolveTranspose (_ksp);
         LIBMESH_CHKERRABORT(ierr);

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         LIBMESH_CHKERRABORT(ierr);

// 2.2.1 & newer style
#else

  Mat submat = NULL;
  Mat subprecond = NULL;
  Vec subrhs = NULL;
  Vec subsolution = NULL;
  VecScatter scatter = NULL;
  PetscMatrix<Number>* subprecond_matrix = NULL;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERRABORT(ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERRABORT(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
	{
	  _create_complement_is(rhs_in);
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;
	  Mat submat1 = NULL;
	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScale(subvec1,-1.0);
	  LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#else
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#endif

	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  LIBMESH_CHKERRABORT(ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond,
						      this->comm());
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

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
         LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolveTranspose (_ksp, subrhs, subsolution);
      LIBMESH_CHKERRABORT(ierr);
    }
  else
    {
      ierr = KSPSolveTranspose (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERRABORT(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         LIBMESH_CHKERRABORT(ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  /* Before we delete subprecond_matrix, we should give the
	     _preconditioner a different matrix.  */
	  this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
	  delete subprecond_matrix;
	  subprecond_matrix = NULL;
	}

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERRABORT(ierr);
    }

#endif

  STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}


template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (const ShellMatrix<T>& shell_matrix,
			     NumericVector<T>& solution_in,
			     NumericVector<T>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{

#if PETSC_VERSION_LESS_THAN(2,3,1)
  // FIXME[JWP]: There will be a bunch of unused variable warnings
  // for older PETScs here.
  libMesh::out << "This method has been developed with PETSc 2.3.1.  "
	        << "No one has made it backwards compatible with older "
	        << "versions of PETSc so far; however, it might work "
	        << "without any change with some older version." << std::endl;
  libmesh_error();
  return std::make_pair(0,0.0);

#else

#if PETSC_VERSION_LESS_THAN(3,1,0)
  if(_restrict_solve_to_is!=NULL)
    {
      libMesh::out << "The current implementation of subset solves with "
		   << "shell matrices requires PETSc version 3.1 or above.  "
		   << "Older PETSc version do not support automatic "
		   << "submatrix generation of shell matrices."
		   << std::endl;
      libmesh_error();
    }
#endif

  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

  this->init ();

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  Mat submat = NULL;
  Vec subrhs = NULL;
  Vec subsolution = NULL;
  VecScatter scatter = NULL;

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
			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
			&mat);
  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERRABORT(ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
	{
	  _create_complement_is(rhs_in);
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;
	  Mat submat1 = NULL;
	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScale(subvec1,-1.0);
	  LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  /* This point can't be reached, see above.  */
	  libmesh_assert(false);
#else
	  ierr = MatGetSubMatrix(mat,
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#endif

	  // The following lines would be correct, but don't work
	  // correctly in PETSc up to 3.1.0-p5.  See discussion in
	  // petsc-users of Nov 9, 2010.
	  //
	  // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  // LIBMESH_CHKERRABORT(ierr);
	  //
	  // We workaround by using a temporary vector.  Note that the
	  // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
	  // so this is no effective performance loss.
	  Vec subvec2 = NULL;
	  ierr = VecCreate(this->comm().get(),&subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = MatMult(submat1,subvec1,subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecAXPY(subrhs,1.0,subvec2);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  LIBMESH_CHKERRABORT(ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, submat,
			     DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);
    }
  else
    {
      ierr = KSPSetOperators(_ksp, mat, mat,
			     DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
  LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERRABORT(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERRABORT(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERRABORT(ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERRABORT(ierr);

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERRABORT(ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  LIBMESH_CHKERRABORT(ierr);

  STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);

#endif

}



template <typename T>
std::pair<unsigned int, Real>
PetscLinearSolver<T>::solve (const ShellMatrix<T>& shell_matrix,
			     const SparseMatrix<T>& precond_matrix,
			     NumericVector<T> &solution_in,
			     NumericVector<T> &rhs_in,
			     const double tol,
			     const unsigned int m_its)
{

#if PETSC_VERSION_LESS_THAN(2,3,1)
  // FIXME[JWP]: There will be a bunch of unused variable warnings
  // for older PETScs here.
  libMesh::out << "This method has been developed with PETSc 2.3.1.  "
	        << "No one has made it backwards compatible with older "
	        << "versions of PETSc so far; however, it might work "
	        << "without any change with some older version." << std::endl;
  libmesh_error();
  return std::make_pair(0,0.0);

#else

#if PETSC_VERSION_LESS_THAN(3,1,0)
  if(_restrict_solve_to_is!=NULL)
    {
      libMesh::out << "The current implementation of subset solves with "
		   << "shell matrices requires PETSc version 3.1 or above.  "
		   << "Older PETSc version do not support automatic "
		   << "submatrix generation of shell matrices."
		   << std::endl;
      libmesh_error();
    }
#endif

  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  const PetscMatrix<T>* precond  = libmesh_cast_ptr<const PetscMatrix<T>*>(&precond_matrix);
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

  this->init ();

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  Mat submat = NULL;
  Mat subprecond = NULL;
  Vec subrhs = NULL;
  Vec subsolution = NULL;
  VecScatter scatter = NULL;
  PetscMatrix<Number>* subprecond_matrix = NULL;

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
			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
			&mat);
  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  LIBMESH_CHKERRABORT(ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERRABORT(ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = MatGetSubMatrix(const_cast<PetscMatrix<T>*>(precond)->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      LIBMESH_CHKERRABORT(ierr);
#endif

      /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
	{
	  _create_complement_is(rhs_in);
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;
	  Mat submat1 = NULL;
	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  LIBMESH_CHKERRABORT(ierr);

	  ierr = VecScale(subvec1,-1.0);
	  LIBMESH_CHKERRABORT(ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  /* This point can't be reached, see above.  */
	  libmesh_assert(false);
#else
	  ierr = MatGetSubMatrix(mat,
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  LIBMESH_CHKERRABORT(ierr);
#endif

	  // The following lines would be correct, but don't work
	  // correctly in PETSc up to 3.1.0-p5.  See discussion in
	  // petsc-users of Nov 9, 2010.
	  //
	  // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  // LIBMESH_CHKERRABORT(ierr);
	  //
	  // We workaround by using a temporary vector.  Note that the
	  // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
	  // so this is no effective performance loss.
	  Vec subvec2 = NULL;
	  ierr = VecCreate(this->comm().get(),&subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecSetFromOptions(subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = MatMult(submat1,subvec1,subvec2);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecAXPY(subrhs,1.0,subvec2);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  LIBMESH_CHKERRABORT(ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond,
						      this->comm());
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix<T>*>(precond)->mat(),
			     DIFFERENT_NONZERO_PATTERN);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
      {
	this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number>&>(precond_matrix));
        this->_preconditioner->init();
      }
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
  LIBMESH_CHKERRABORT(ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      LIBMESH_CHKERRABORT(ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      LIBMESH_CHKERRABORT(ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  LIBMESH_CHKERRABORT(ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  LIBMESH_CHKERRABORT(ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  LIBMESH_CHKERRABORT(ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      LIBMESH_CHKERRABORT(ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      LIBMESH_CHKERRABORT(ierr);

      if(this->_preconditioner)
	{
	  /* Before we delete subprecond_matrix, we should give the
	     _preconditioner a different matrix.  */
	  this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number>&>(precond_matrix));
          this->_preconditioner->init();
	  delete subprecond_matrix;
	  subprecond_matrix = NULL;
	}

      ierr = LibMeshVecDestroy(&subsolution);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&submat);
      LIBMESH_CHKERRABORT(ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      LIBMESH_CHKERRABORT(ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  LIBMESH_CHKERRABORT(ierr);

  STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);

#endif

}



template <typename T>
void PetscLinearSolver<T>::get_residual_history(std::vector<double>& hist)
{
  PetscErrorCode ierr = 0;
  PetscInt its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  LIBMESH_CHKERRABORT(ierr);

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
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  LIBMESH_CHKERRABORT(ierr);

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
      ierr = KSPSetType (_ksp, (char*) KSPCG);         LIBMESH_CHKERRABORT(ierr); return;

    case CR:
      ierr = KSPSetType (_ksp, (char*) KSPCR);         LIBMESH_CHKERRABORT(ierr); return;

    case CGS:
      ierr = KSPSetType (_ksp, (char*) KSPCGS);        LIBMESH_CHKERRABORT(ierr); return;

    case BICG:
      ierr = KSPSetType (_ksp, (char*) KSPBICG);       LIBMESH_CHKERRABORT(ierr); return;

    case TCQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTCQMR);      LIBMESH_CHKERRABORT(ierr); return;

    case TFQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTFQMR);      LIBMESH_CHKERRABORT(ierr); return;

    case LSQR:
      ierr = KSPSetType (_ksp, (char*) KSPLSQR);       LIBMESH_CHKERRABORT(ierr); return;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, (char*) KSPBCGS);       LIBMESH_CHKERRABORT(ierr); return;

    case MINRES:
      ierr = KSPSetType (_ksp, (char*) KSPMINRES);     LIBMESH_CHKERRABORT(ierr); return;

    case GMRES:
      ierr = KSPSetType (_ksp, (char*) KSPGMRES);      LIBMESH_CHKERRABORT(ierr); return;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, (char*) KSPRICHARDSON); LIBMESH_CHKERRABORT(ierr); return;

    case CHEBYSHEV:
#if defined(LIBMESH_HAVE_PETSC) && PETSC_VERSION_LESS_THAN(3,3,0)
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV);  LIBMESH_CHKERRABORT(ierr); return;
#else
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYSHEV);  LIBMESH_CHKERRABORT(ierr); return;
#endif


    default:
      libMesh::err << "ERROR:  Unsupported PETSC Solver: "
		    << Utility::enum_to_string(this->_solver_type) << std::endl
		    << "Continuing with PETSC defaults" << std::endl;
    }
}

template <typename T>
void PetscLinearSolver<T>::print_converged_reason()
{
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);
  libMesh::out << "Linear solver convergence/divergence reason: " << KSPConvergedReasons[reason] << std::endl;
}



template <typename T>
PetscErrorCode PetscLinearSolver<T>::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
{
  /* Get the matrix context.  */
  PetscErrorCode ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
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
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
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
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
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
