// $Id: petsc_interface.C,v 1.12 2003-03-04 22:31:16 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#include "mesh_common.h"

#ifdef HAVE_PETSC


// C++ includes

// Local Includes
#include "petsc_interface.h"



/*----------------------- functions ----------------------------------*/
template <typename T>
void PetscInterface<T>::clear ()
{
  if (initialized())
    {
      _is_initialized = false;

      int ierr=0;

      ierr = SLESDestroy(_sles); CHKERRQ(ierr);

      // Mimic PETSc default solver and preconditioner
      _solver_type         = GMRES;

      if (libMeshBase::n_processors() == 1)
	_preconditioner_type = ILU_PRECOND;
      else
	_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
}



template <typename T>
void PetscInterface<T>::init ()
{
  int ierr=0;
  
  // Initialize the data structures if not done so already.
  if (!initialized())
    {
      _is_initialized = true;

      // Create the linear solver context
      ierr = SLESCreate (PETSC_COMM_WORLD, &_sles); CHKERRQ(ierr);
      
      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (_sles, &_ksp); CHKERRQ(ierr);
      ierr = SLESGetPC        (_sles, &_pc);  CHKERRQ(ierr);
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      // 2.1.0 - style
      //ierr = KSPSetInitialGuessNonzero (ksp); CHKERRQ(ierr);
      // 2.1.1 - style
      ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE); CHKERRQ(ierr);
      
      // Set user-specified  solver and preconditioner types
      set_petsc_solver_type();
      set_petsc_preconditioner_type();
      
      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.
      
      ierr = SLESSetFromOptions (_sles);                   CHKERRQ(ierr);
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
  init ();
  
  PetscMatrix<T>& matrix   = dynamic_cast<PetscMatrix<T>&>(matrix_in);
  PetscVector<T>& solution = dynamic_cast<PetscVector<T>&>(solution_in);
  PetscVector<T>& rhs      = dynamic_cast<PetscVector<T>&>(rhs_in);
  
  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;

  // Close the matrix and vectors in case this wasn't already done.
  matrix.close ();
  solution.close ();
  rhs.close ();

  
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = SLESSetOperators(_sles, matrix.mat, matrix.mat,
			  SAME_NONZERO_PATTERN);             CHKERRQ(ierr);

  
  ierr = KSPSetTolerances (_ksp, tol, tol,
			   PETSC_DEFAULT,
			   max_its);                         CHKERRQ(ierr);

  
  // Solve the linear system
  ierr = SLESSolve (_sles, rhs.vec, solution.vec, &its);     CHKERRQ(ierr);
  
  
  ierr = KSPGetResidualNorm (_ksp, &final_resid);            CHKERRQ(ierr);
  
  std::pair<unsigned int, Real> p (its, final_resid);

  return p;
}



template <typename T>
void PetscInterface<T>::set_petsc_solver_type()
{
  int ierr = 0;
  
  switch (_solver_type)
    {

    case CG:
      ierr = KSPSetType (_ksp, (char*) KSPCG); CHKERRQ(ierr); return;

    case CR:
      ierr = KSPSetType (_ksp, (char*) KSPCR); CHKERRQ(ierr); return;

    case CGS:
      ierr = KSPSetType (_ksp, (char*) KSPCGS); CHKERRQ(ierr); return;

    case BICG:
      ierr = KSPSetType (_ksp, (char*) KSPBICG); CHKERRQ(ierr); return;

    case TCQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTCQMR); CHKERRQ(ierr); return;

    case TFQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTFQMR); CHKERRQ(ierr); return;

    case LSQR:
      ierr = KSPSetType (_ksp, (char*) KSPLSQR); CHKERRQ(ierr); return;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, (char*) KSPBCGS); CHKERRQ(ierr); return;

    case MINRES:
      ierr = KSPSetType (_ksp, (char*) KSPMINRES); CHKERRQ(ierr); return;

    case GMRES:
      ierr = KSPSetType (_ksp, (char*) KSPGMRES); CHKERRQ(ierr); return;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, (char*) KSPRICHARDSON); CHKERRQ(ierr); return;

    case CHEBYSHEV:
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV); CHKERRQ(ierr); return;

    default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
		<< _solver_type                     << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
}



template <typename T>
void PetscInterface<T>::set_petsc_preconditioner_type()
{
  int ierr = 0;
 
  switch (_preconditioner_type)
    {
    case IDENTITY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCNONE); CHKERRQ(ierr); return;
	
    case CHOLESKY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCCHOLESKY); CHKERRQ(ierr); return;

    case ICC_PRECOND:
      ierr = PCSetType (_pc, (char*) PCICC); CHKERRQ(ierr); return;

    case ILU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCILU); CHKERRQ(ierr); return;

    case LU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCLU); CHKERRQ(ierr); return;
      
    case ASM_PRECOND:
      ierr = PCSetType (_pc, (char*) PCASM); CHKERRQ(ierr); return;

    case JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCJACOBI); CHKERRQ(ierr); return;

    case BLOCK_JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCBJACOBI); CHKERRQ(ierr); return;

    case SOR_PRECOND:
      ierr = PCSetType (_pc, (char*) PCSOR); CHKERRQ(ierr); return;

    case EISENSTAT_PRECOND:
      ierr = PCSetType (_pc, (char*) PCEISENSTAT); CHKERRQ(ierr); return;

    default:
      std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
		<< _preconditioner_type             << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscInterface<Number>;
 


#endif // #ifdef HAVE_PETSC
