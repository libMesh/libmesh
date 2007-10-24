// $Id: petsc_interface.C,v 1.24 2004-04-17 03:02:50 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

#ifdef HAVE_PETSC


// C++ includes

// Local Includes
#include "petsc_interface.h"



/*----------------------- functions ----------------------------------*/
template <typename T>
void PetscInterface<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

// 2.1.x & earlier style      
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
      
      ierr = SLESDestroy(_sles);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);

// 2.2.0 & newer style
#else
      
      ierr = KSPDestroy(_ksp);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
#endif
	     
      // Mimic PETSc default solver and preconditioner
      this->_solver_type           = GMRES;

      if (libMesh::n_processors() == 1)
	this->_preconditioner_type = ILU_PRECOND;
      else
	this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
}



template <typename T>
void PetscInterface<T>::init ()
{
  int ierr=0;
  
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

// 2.1.x & earlier style      
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
      
      // Create the linear solver context
      ierr = SLESCreate (PETSC_COMM_WORLD, &_sles);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      ierr = SLESGetPC        (_sles, &_pc);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();
      this->set_petsc_preconditioner_type();
      
      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.
      
      ierr = SLESSetFromOptions (_sles);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);

// 2.2.0 & newer style
#else
      
      // Create the linear solver context
      ierr = KSPCreate (PETSC_COMM_WORLD, &_ksp);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Create the preconditioner context
      ierr = KSPGetPC        (_ksp, &_pc);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type();
      this->set_petsc_preconditioner_type();
      
      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.
      
      ierr = KSPSetFromOptions (_ksp);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);

	       
#endif
	       
    }
}



template <typename T>
std::pair<unsigned int, Real> 
PetscInterface<T>::solve (SparseMatrix<T> &matrix_in,
			  NumericVector<T> &solution_in,
			  NumericVector<T> &rhs_in,
			  const double tol,
			  const unsigned int m_its)
{
  this->init ();
  
  PetscMatrix<T>* matrix   = dynamic_cast<PetscMatrix<T>*>(&matrix_in);
  PetscVector<T>* solution = dynamic_cast<PetscVector<T>*>(&solution_in);
  PetscVector<T>* rhs      = dynamic_cast<PetscVector<T>*>(&rhs_in);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(matrix   != NULL);
  assert(solution != NULL);
  assert(rhs      != NULL);
  
  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;

  // Close the matrix and vectors in case this wasn't already done.
  matrix->close ();
  solution->close ();
  rhs->close ();

  
// 2.1.x & earlier style      
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
      
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = SLESSetOperators(_sles, matrix->mat, matrix->mat,
			  SAME_NONZERO_PATTERN);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);


  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);


  // Solve the linear system
  ierr = SLESSolve (_sles, rhs->vec, solution->vec, &its);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);


  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

// 2.2.0 & newer style
#else
      
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp, matrix->mat, matrix->mat,
			 SAME_NONZERO_PATTERN);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);


  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
   
  // Solve the linear system
  ierr = KSPSolve (_ksp);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
#endif

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}



template <typename T>
void PetscInterface<T>::set_petsc_solver_type()
{
  int ierr = 0;
  
  switch (this->_solver_type)
    {

    case CG:
      ierr = KSPSetType (_ksp, (char*) KSPCG);         CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case CR:
      ierr = KSPSetType (_ksp, (char*) KSPCR);         CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case CGS:
      ierr = KSPSetType (_ksp, (char*) KSPCGS);        CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case BICG:
      ierr = KSPSetType (_ksp, (char*) KSPBICG);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case TCQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTCQMR);      CHKERRABORT(PETSC_COMM_WORLD,ierr); return;
 
    case TFQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTFQMR);      CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case LSQR:
      ierr = KSPSetType (_ksp, (char*) KSPLSQR);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, (char*) KSPBCGS);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case MINRES:
      ierr = KSPSetType (_ksp, (char*) KSPMINRES);     CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case GMRES:
      ierr = KSPSetType (_ksp, (char*) KSPGMRES);      CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, (char*) KSPRICHARDSON); CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case CHEBYSHEV: 
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV);  CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
		<< this->_solver_type               << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
}



template <typename T>
void PetscInterface<T>::set_petsc_preconditioner_type()
{
  int ierr = 0;
 
  switch (this->_preconditioner_type)
    {
    case IDENTITY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCNONE);      CHKERRABORT(PETSC_COMM_WORLD,ierr); return;
	
    case CHOLESKY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCCHOLESKY);  CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case ICC_PRECOND:
      ierr = PCSetType (_pc, (char*) PCICC);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case ILU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCILU);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case LU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCLU);        CHKERRABORT(PETSC_COMM_WORLD,ierr); return;
      
    case ASM_PRECOND:
      ierr = PCSetType (_pc, (char*) PCASM);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCJACOBI);    CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case BLOCK_JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCBJACOBI);   CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case SOR_PRECOND:
      ierr = PCSetType (_pc, (char*) PCSOR);       CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    case EISENSTAT_PRECOND:
      ierr = PCSetType (_pc, (char*) PCEISENSTAT); CHKERRABORT(PETSC_COMM_WORLD,ierr); return;

    default:
      std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
		<< this->_preconditioner_type       << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscInterface<Number>;
 


#endif // #ifdef HAVE_PETSC
