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

#include "libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC


// C++ includes

// Local Includes
#include "libmesh_logging.h"
#include "petsc_linear_solver.h"
#include "shell_matrix.h"
#include "petsc_preconditioner.h"

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
    preconditioner->init();

    return 0;
  }
  

  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)
  {
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

    PetscVector<Number> x_vec(x);
    PetscVector<Number> y_vec(y);

    preconditioner->apply(x_vec,y_vec);

    return 0;
  }
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC pc)
  {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
    preconditioner->init();

    return 0;
  }

  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y)
  {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);

    PetscVector<Number> x_vec(x);
    PetscVector<Number> y_vec(y);

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
      this->_is_initialized = false;

      int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,2,0)

  // 2.1.x & earlier style      
  ierr = SLESDestroy(_sles);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else

  // 2.2.0 & newer style
  ierr = KSPDestroy(_ksp);
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
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
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

      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
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
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
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
      ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),SAME_NONZERO_PATTERN);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);       

      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
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
        PCShellSetContext(_pc,(void*)this->_preconditioner);
        PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
        PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
      }
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

  if(this->_preconditioner)
    this->_preconditioner->set_matrix(matrix_in);

  
// 2.1.x & earlier style      
#if PETSC_VERSION_LESS_THAN(2,2,0)
      
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = SLESSetOperators(_sles, matrix->mat(), precond->mat(),
			  SAME_NONZERO_PATTERN);
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
      
  // Set operators. The input matrix works as the preconditioning matrix
  if(!this->same_preconditioner)
  {
    ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
	  		   SAME_NONZERO_PATTERN);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  }
  else
  {
    ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
               SAME_PRECONDITIONER);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
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
  std::cout << "This method has been developed with PETSc 2.3.1.  "
	    << "No one has made it backwards compatible with older "
	    << "versions of PETSc so far; however, it might work "
	    << "without any change with some older version." << std::endl;
  libmesh_error();
  return std::make_pair(0,0.0);

#else

  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

  this->init ();

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;

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
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp, mat, mat,
			 SAME_NONZERO_PATTERN);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Destroy the matrix.
  ierr = MatDestroy(mat);
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
  std::cout << "This method has been developed with PETSc 2.3.1.  "
	    << "No one has made it backwards compatible with older "
	    << "versions of PETSc so far; however, it might work "
	    << "without any change with some older version." << std::endl;
  libmesh_error();
  return std::make_pair(0,0.0);

#else

  START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  const PetscMatrix<T>* precond  = libmesh_cast_ptr<const PetscMatrix<T>*>(&precond_matrix);
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = libmesh_cast_ptr<PetscVector<T>*>(&rhs_in);

  this->init ();

  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;

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
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix<T>*>(precond)->mat(),
			 DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  if(this->_preconditioner)
    this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number>&>(precond_matrix));

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Destroy the matrix.
  ierr = MatDestroy(mat);
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
      std::cerr << "No iterations have been performed, returning 0." << std::endl;
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
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
		<< this->_solver_type               << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
}

template <typename T>
void PetscLinearSolver<T>::print_converged_reason()
{
#if PETSC_VERSION_LESS_THAN(2,3,1)
  std::cout << "This method is currently not supported "
	    << "(but may work!) for Petsc 2.3.0 and earlier." << std::endl;
#else
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);
  
  //  KSP_CONVERGED_RTOL (residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side)
  //  KSP_CONVERGED_ATOL (residual 2-norm less than abstol)
  //  KSP_CONVERGED_ITS (used by the preonly preconditioner that always uses ONE iteration) 
  //  KSP_CONVERGED_STEP_LENGTH
  //  KSP_DIVERGED_ITS  (required more than its to reach convergence)
  //  KSP_DIVERGED_DTOL (residual norm increased by a factor of divtol)
  //  KSP_DIVERGED_NAN (residual norm became Not-a-number likely do to 0/0)
  //  KSP_DIVERGED_BREAKDOWN (generic breakdown in method)

  switch (reason)
    {
    case KSP_CONVERGED_RTOL:
       {
	std::cout << "Linear solver converged, relative tolerance reached." << std::endl;
	break;
       }
    case KSP_CONVERGED_ATOL:
       {
	 std::cout << "Linear solver converged, absolute tolerance reached." << std::endl;
	 break;
       }

      // Divergence
    case KSP_DIVERGED_ITS:
       {
	 std::cout << "Linear solver diverged, max no. of iterations reached." << std::endl;
	 break;
       }
    case KSP_DIVERGED_DTOL:
       {
	 std::cout << "Linear solver diverged, residual norm increase by dtol (default 1.e5)." << std::endl;
	 break;
       }
    case KSP_DIVERGED_NAN:
       {
	 std::cout << "Linear solver diverged, residual norm is NaN." << std::endl;
	 break;
       }
    case KSP_DIVERGED_BREAKDOWN:
       {
	 std::cout << "Linear solver diverged, generic breakdown in the method." << std::endl;
	 break;
       }
    default:
      {
	std::cout << "Unknown/unsupported con(di)vergence reason: " << reason << std::endl;
      }
    }
#endif
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
  PetscVector<T> arg_global(arg);
  PetscVector<T> dest_global(dest);

  /* Call the user function.  */
  shell_matrix.vector_mult(dest_global,arg_global);

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
  PetscVector<T> dest_global(dest);

  /* Call the user function.  */
  shell_matrix.get_diagonal(dest_global);

  return ierr;
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscLinearSolver<Number>;
 


#endif // #ifdef LIBMESH_HAVE_PETSC
