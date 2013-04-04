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

namespace libMesh
{

extern "C"
{
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif


#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
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

    PetscVector<Number> x_vec(x, preconditioner->communicator());
    PetscVector<Number> y_vec(y, preconditioner->communicator());

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

    PetscVector<Number> x_vec(x/*, preconditioner->communicator()*/);
    PetscVector<Number> y_vec(y/*, preconditioner->communicator()*/);

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
	  int ierr = LibMeshISDestroy(&_restrict_solve_to_is);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  _restrict_solve_to_is = NULL;
	}

      if(_restrict_solve_to_is_complement!=NULL)
	{
	  int ierr = LibMeshISDestroy(&_restrict_solve_to_is_complement);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  _restrict_solve_to_is_complement = NULL;
	}

      this->_is_initialized = false;

      int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,2,0)

  // 2.1.x & earlier style
  ierr = SLESDestroy(_sles);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else

  // 2.2.0 & newer style
  ierr = LibMeshKSPDestroy(&_ksp);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif


      // Mimic PETSc default solver and preconditioner
      this->_solver_type           = GMRES;

      if(!this->_preconditioner)
      {
        if (libMesh::n_processors() == 1)
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

      int ierr=0;

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

      // Create the linear solver context
      ierr = SLESCreate (libMesh::COMM_WORLD, &_sles);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = SLESGetPC        (_sles, &_pc);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = SLESSetFromOptions (_sles);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.0 & newer style
#else

      // Create the linear solver context
      ierr = KSPCreate (libMesh::COMM_WORLD, &_ksp);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //CHKERRABORT(libMesh::COMM_WORLD,ierr);


#endif

      // Have the Krylov subspace method use our good initial guess
      // rather than 0, unless the user requested a KSPType of
      // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_VERSION_RELEASE
      // Pre-3.0 and petsc-dev (as of October 2012) use non-const versions
      KSPType ksp_type;
#else
      const KSPType ksp_type;
#endif

      ierr = KSPGetType (_ksp, &ksp_type);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
                 CHKERRABORT(libMesh::COMM_WORLD,ierr);
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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

      int ierr=0;

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

      // Create the linear solver context
      ierr = SLESCreate (libMesh::COMM_WORLD, &_sles);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = SLESGetPC        (_sles, &_pc);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = SLESSetFromOptions (_sles);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.0 & newer style
#else

      // Create the linear solver context
      ierr = KSPCreate (libMesh::COMM_WORLD, &_ksp);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);


      //ierr = PCCreate (libMesh::COMM_WORLD, &_pc);
      //     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set operators. The input matrix works as the preconditioning matrix
      ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),DIFFERENT_NONZERO_PATTERN);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();

      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.

      ierr = KSPSetFromOptions (_ksp);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      // NOT NECESSARY!!!!
      //ierr = PCSetFromOptions (_pc);
      //CHKERRABORT(libMesh::COMM_WORLD,ierr);


#endif

      // Have the Krylov subspace method use our good initial guess
      // rather than 0, unless the user requested a KSPType of
      // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_VERSION_LESS_THAN(3,4,0) || !PETSC_VERSION_RELEASE
      KSPType ksp_type;
#else
      const KSPType ksp_type;
#endif

      ierr = KSPGetType (_ksp, &ksp_type);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      if (strcmp(ksp_type, "preonly"))
        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
                 CHKERRABORT(libMesh::COMM_WORLD,ierr);
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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  int ierr=0;

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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for(size_t i=0; i<dofs->size(); i++)
	{
	  petsc_dofs[i] = (*dofs)[i];
	}

      ierr = ISCreateLibMesh(libMesh::COMM_WORLD,dofs->size(),petsc_dofs,PETSC_OWN_POINTER,&_restrict_solve_to_is);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
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
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


  // Solve the linear system
  ierr = SLESSolve (_sles, rhs->vec(), solution->vec(), &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.0
#elif PETSC_VERSION_LESS_THAN(2,2,1)

  if(_restrict_solve_to_is!=NULL)
    {
      libmesh_not_implemented();
    }

  // Set operators. The input matrix works as the preconditioning matrix
  //ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
	 //			 SAME_NONZERO_PATTERN);
	 //       CHKERRABORT(libMesh::COMM_WORLD,ierr);


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
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolve (_ksp);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

      ierr = VecCreate(libMesh::COMM_WORLD,&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecCreate(libMesh::COMM_WORLD,&subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond);
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
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
  ierr = SLESCreate (libMesh::COMM_WORLD, &sles);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SLESSetOperators (sles, matrix->mat(), precond->mat(), this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  KSP ksp;
  ierr = SLESGetKSP (sles, &ksp);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SLESSetUp (sles, rhs->vec(), solution->vec());
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // See http://tccc.iesl.forth.gr/AMS_EPEAEK/Elements/doc/in_html/petsc/KSP/KSPSolveTrans.html#KSPSolveTrans
  ierr = SLESSolveTrans (ksp, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  CHKERRABORT(libMesh::COMM_WORLD,ierr);


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
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolveTranspose (_ksp);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

      ierr = VecCreate(libMesh::COMM_WORLD,&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecCreate(libMesh::COMM_WORLD,&subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond);
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			     this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolveTranspose (_ksp, subrhs, subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      ierr = KSPSolveTranspose (_ksp, rhs->vec(), solution->vec());
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
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
  ierr = MatCreateShell(libMesh::COMM_WORLD,
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

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(libMesh::COMM_WORLD,&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecCreate(libMesh::COMM_WORLD,&subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  /* This point can't be reached, see above.  */
	  libmesh_assert(false);
#else
	  ierr = MatGetSubMatrix(mat,
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

	  // The following lines would be correct, but don't work
	  // correctly in PETSc up to 3.1.0-p5.  See discussion in
	  // petsc-users of Nov 9, 2010.
	  //
	  // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  // CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  //
	  // We workaround by using a temporary vector.  Note that the
	  // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
	  // so this is no effective performance loss.
	  Vec subvec2 = NULL;
	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = MatMult(submat1,subvec1,subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecAXPY(subrhs,1.0,subvec2);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, submat,
			     DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      ierr = KSPSetOperators(_ksp, mat, mat,
			     DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = LibMeshVecDestroy(&subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
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
  ierr = MatCreateShell(libMesh::COMM_WORLD,
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

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT_ADD,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult_add));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Restrict rhs and solution vectors and set operators.  The input
  // matrix works as the preconditioning matrix.
  if(_restrict_solve_to_is!=NULL)
    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(libMesh::COMM_WORLD,&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecCreate(libMesh::COMM_WORLD,&subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      /* This point can't be reached, see above.  */
      libmesh_assert(false);
#else
      ierr = MatGetSubMatrix(mat,
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(const_cast<PetscMatrix<T>*>(precond)->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  /* This point can't be reached, see above.  */
	  libmesh_assert(false);
#else
	  ierr = MatGetSubMatrix(mat,
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

	  // The following lines would be correct, but don't work
	  // correctly in PETSc up to 3.1.0-p5.  See discussion in
	  // petsc-users of Nov 9, 2010.
	  //
	  // ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);
	  // CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  //
	  // We workaround by using a temporary vector.  Note that the
	  // fix in PETsc 3.1.0-p6 uses a temporary vector internally,
	  // so this is no effective performance loss.
	  Vec subvec2 = NULL;
	  ierr = VecCreate(libMesh::COMM_WORLD,&subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec2,is_local_size,PETSC_DECIDE);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = MatMult(submat1,subvec1,subvec2);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecAXPY(subrhs,1.0,subvec2);

	  ierr = LibMeshVecScatterDestroy(&scatter1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshVecDestroy(&subvec1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = LibMeshMatDestroy(&submat1);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond,
			     DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      if(this->_preconditioner)
	{
	  subprecond_matrix = new PetscMatrix<Number>(subprecond);
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else
    {
      ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix<T>*>(precond)->mat(),
			     DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL)
    {
      ierr = KSPSolve (_ksp, subrhs, subsolution);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  if(_restrict_solve_to_is!=NULL)
    {
      switch(_subset_solve_mode)
	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = LibMeshVecScatterDestroy(&scatter);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshVecDestroy(&subrhs);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&submat);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = LibMeshMatDestroy(&subprecond);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Destroy the matrix.
  ierr = LibMeshMatDestroy(&mat);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);

#endif

}



template <typename T>
void PetscLinearSolver<T>::get_residual_history(std::vector<double>& hist)
{
  int ierr = 0;
  int its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Check for early return
  if (its == 0) return;

  // Create space to store the result
  hist.resize(its);

  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i)
    {
      hist[i] = *p;
      p++;
    }
}




template <typename T>
Real PetscLinearSolver<T>::get_initial_residual()
{
  int ierr = 0;
  int its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  int ierr = 0;

  switch (this->_solver_type)
    {

    case CG:
      ierr = KSPSetType (_ksp, (char*) KSPCG);         CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case CR:
      ierr = KSPSetType (_ksp, (char*) KSPCR);         CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case CGS:
      ierr = KSPSetType (_ksp, (char*) KSPCGS);        CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case BICG:
      ierr = KSPSetType (_ksp, (char*) KSPBICG);       CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case TCQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTCQMR);      CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case TFQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTFQMR);      CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case LSQR:
      ierr = KSPSetType (_ksp, (char*) KSPLSQR);       CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, (char*) KSPBCGS);       CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case MINRES:
      ierr = KSPSetType (_ksp, (char*) KSPMINRES);     CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case GMRES:
      ierr = KSPSetType (_ksp, (char*) KSPGMRES);      CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, (char*) KSPRICHARDSON); CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    case CHEBYSHEV:
#if defined(LIBMESH_HAVE_PETSC) && PETSC_VERSION_LESS_THAN(3,3,0)
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
#else
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYSHEV);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
#endif


    default:
      libMesh::err << "ERROR:  Unsupported PETSC Solver: "
		    << this->_solver_type               << std::endl
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
  int ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);

  /* Make \p NumericVector instances around the vectors.  */
  PetscVector<T> arg_global(arg/*, shell_matrix.communicator()*/);
  PetscVector<T> dest_global(dest/*, shell_matrix.communicator()*/);

  /* Call the user function.  */
  shell_matrix.vector_mult(dest_global,arg_global);

  return ierr;
}



template <typename T>
PetscErrorCode PetscLinearSolver<T>::_petsc_shell_matrix_mult_add(Mat mat, Vec arg, Vec add, Vec dest)
{
  /* Get the matrix context.  */
  int ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);

  /* Make \p NumericVector instances around the vectors.  */
  PetscVector<T> arg_global(arg/*, shell_matrix.communicator()*/);
  PetscVector<T> dest_global(dest/*, shell_matrix.communicator()*/);
  PetscVector<T> add_global(add/*, shell_matrix.communicator()*/);

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
  int ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);

  /* Make \p NumericVector instances around the vector.  */
  PetscVector<T> dest_global(dest/*, shell_matrix.communicator()*/);

  /* Call the user function.  */
  shell_matrix.get_diagonal(dest_global);

  return ierr;
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscLinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
