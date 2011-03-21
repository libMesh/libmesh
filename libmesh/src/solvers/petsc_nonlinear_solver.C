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
#include "nonlinear_implicit_system.h"
#include "petsc_nonlinear_solver.h"
#include "petsc_linear_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "dof_map.h"
#include "preconditioner.h"

namespace libMesh
{

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
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
  
  //-------------------------------------------------------------------
  // this function is called by PETSc at the end of each nonlinear step  
  PetscErrorCode
  __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *)
  {
    //int ierr=0;

    //if (its > 0)
      libMesh::out << "  NL step "
		    << std::setw(2) << its
		    << std::scientific
		    << ", |residual|_2 = " << fnorm
		    << std::endl;

    //return ierr;
    return 0;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __libmesh_petsc_snes_residual (SNES, Vec x, Vec r, void *ctx)
  {
    START_LOG("residual()", "PetscNonlinearSolver");
    int ierr=0;

    libmesh_assert (x   != NULL);
    libmesh_assert (r   != NULL);
    libmesh_assert (ctx != NULL);
    
    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);
    
    NonlinearImplicitSystem &sys = solver->system();

    PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
    PetscVector<Number>& R_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.rhs);
    PetscVector<Number> X_global(x), R(r);

    // Use the systems update() to get a good local version of the parallel solution
    X_global.swap(X_sys);
    R.swap(R_sys);

    sys.get_dof_map().enforce_constraints_exactly(sys);
    
    sys.update();

    //Swap back
    X_global.swap(X_sys);
    R.swap(R_sys);

    R.zero();

    if      (solver->residual != NULL) solver->residual (*sys.current_local_solution.get(), R, sys);
    else if (solver->matvec   != NULL) solver->matvec   (*sys.current_local_solution.get(), &R, NULL, sys);
    else libmesh_error();
    
    R.close();
    X_global.close();

    STOP_LOG("residual()", "PetscNonlinearSolver");
    
    return ierr;
  }


  
  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __libmesh_petsc_snes_jacobian (SNES, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    START_LOG("jacobian()", "PetscNonlinearSolver");
    int ierr=0;
    
    libmesh_assert (ctx != NULL);
    
    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);

    NonlinearImplicitSystem &sys = solver->system();
    
    PetscMatrix<Number> PC(*pc);
    PetscMatrix<Number> Jac(*jac);
    PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
    PetscMatrix<Number>& Jac_sys = *libmesh_cast_ptr<PetscMatrix<Number>*>(sys.matrix);
    PetscVector<Number> X_global(x);

    // Set the dof maps
    PC.attach_dof_map(sys.get_dof_map());
    Jac.attach_dof_map(sys.get_dof_map());

    // Use the systems update() to get a good local version of the parallel solution
    X_global.swap(X_sys);
    Jac.swap(Jac_sys);

    sys.get_dof_map().enforce_constraints_exactly(sys);
    sys.update();

    X_global.swap(X_sys);
    Jac.swap(Jac_sys);

    PC.zero();

    if      (solver->jacobian != NULL) solver->jacobian (*sys.current_local_solution.get(), PC, sys);
    else if (solver->matvec   != NULL) solver->matvec   (*sys.current_local_solution.get(), NULL, &PC, sys);
    else libmesh_error();
    
    PC.close();
    Jac.close();
    X_global.close();
    
    *msflag = SAME_NONZERO_PATTERN;
    
    STOP_LOG("jacobian()", "PetscNonlinearSolver");

    return ierr;
  }
    
} // end extern "C"
//---------------------------------------------------------------------



//---------------------------------------------------------------------
// PetscNonlinearSolver<> methods
template <typename T>
void PetscNonlinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

      ierr = SNESDestroy(_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
}



template <typename T>
void PetscNonlinearSolver<T>::init ()
{  
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;
      
      int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,1,2)
      // At least until Petsc 2.1.1, the SNESCreate had a different calling syntax.
      // The second argument was of type SNESProblemType, and could have a value of
      // either SNES_NONLINEAR_EQUATIONS or SNES_UNCONSTRAINED_MINIMIZATION.
      ierr = SNESCreate(libMesh::COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else

      ierr = SNESCreate(libMesh::COMM_WORLD,&_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif

             
#if PETSC_VERSION_LESS_THAN(2,3,3)
      ierr = SNESSetMonitor (_snes, __libmesh_petsc_snes_monitor,
			     this, PETSC_NULL);
#else
      // API name change in PETSc 2.3.3
      ierr = SNESMonitorSet (_snes, __libmesh_petsc_snes_monitor,
			     this, PETSC_NULL);
#endif
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      // Cannot call SNESSetOptions before SNESSetFunction when using
      // any matrix free options with PETSc 3.1.0+
      ierr = SNESSetFromOptions(_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif 
             
      if(this->_preconditioner)
      {
        KSP ksp;	 
        ierr = SNESGetKSP (_snes, &ksp);
               CHKERRABORT(libMesh::COMM_WORLD,ierr);
        PC pc;
        ierr = KSPGetPC(ksp,&pc);
               CHKERRABORT(libMesh::COMM_WORLD,ierr);
               
        PCSetType(pc, PCSHELL);
        PCShellSetContext(pc,(void*)this->_preconditioner);

        //Re-Use the shell functions from petsc_linear_solver
        PCShellSetSetUp(pc,__libmesh_petsc_preconditioner_setup);
        PCShellSetApply(pc,__libmesh_petsc_preconditioner_apply);
      }
    }
}



template <typename T>
std::pair<unsigned int, Real> 
PetscNonlinearSolver<T>::solve (SparseMatrix<T>&  jac_in,  // System Jacobian Matrix
				NumericVector<T>& x_in,    // Solution vector
				NumericVector<T>& r_in,    // Residual vector
				const double,              // Stopping tolerance
				const unsigned int) 
{
  START_LOG("solve()", "PetscNonlinearSolver");
  this->init ();
  
  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T>* jac = libmesh_cast_ptr<PetscMatrix<T>*>(&jac_in);
  PetscVector<T>* x   = libmesh_cast_ptr<PetscVector<T>*>(&x_in);
  PetscVector<T>* r   = libmesh_cast_ptr<PetscVector<T>*>(&r_in);

  int ierr=0;
  int n_iterations =0;
  // Should actually be a PetscReal, but I don't know which version of PETSc first introduced PetscReal
  Real final_residual_norm=0.; 

  ierr = SNESSetFunction (_snes, r->vec(), __libmesh_petsc_snes_residual, this);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SNESSetJacobian (_snes, jac->mat(), jac->mat(), __libmesh_petsc_snes_jacobian, this);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

   // Have the Krylov subspace method use our good initial guess rather than 0
   KSP ksp;	 
   ierr = SNESGetKSP (_snes, &ksp);
          CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values
  ierr = KSPSetTolerances (ksp, this->initial_linear_tolerance, PETSC_DEFAULT,
                           PETSC_DEFAULT, this->max_linear_iterations);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set the tolerances for the non-linear solver.
  ierr = SNESSetTolerances(_snes, this->absolute_residual_tolerance, this->relative_residual_tolerance,
                           this->absolute_step_tolerance, this->max_nonlinear_iterations, this->max_function_evaluations);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  //Pull in command-line options       
  KSPSetFromOptions(ksp);
  SNESSetFromOptions(_snes);

  //Set the preconditioning matrix
  if(this->_preconditioner)
    this->_preconditioner->set_matrix(jac_in);

//    ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
//           CHKERRABORT(libMesh::COMM_WORLD,ierr);
      	 
// Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
// the last one being a pointer to an int to hold the number of iterations required.
# if PETSC_VERSION_LESS_THAN(2,2,0)

 ierr = SNESSolve (_snes, x->vec(), &n_iterations);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.x style	
#elif PETSC_VERSION_LESS_THAN(2,3,0)
	
 ierr = SNESSolve (_snes, x->vec());
        CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.3.x & newer style	
#else
	
  ierr = SNESSolve (_snes, PETSC_NULL, x->vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SNESGetIterationNumber(_snes,&n_iterations);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
  ierr = SNESGetFunctionNorm(_snes,&final_residual_norm);
	 CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
#endif

  // Get and store the reason for convergence
  SNESGetConvergedReason(_snes, &_reason);

  //Based on Petsc 2.3.3 documentation all diverged reasons are negative
  this->converged = (_reason >= 0);

  this->clear();
  
  STOP_LOG("solve()", "PetscNonlinearSolver");
		 
  // return the # of its. and the final residual norm.
  return std::make_pair(n_iterations, final_residual_norm);
}



template <typename T>
void PetscNonlinearSolver<T>::print_converged_reason()
{
  
  libMesh::out << "Nonlinear solver convergence/divergence reason: "
               << SNESConvergedReasons[this->get_converged_reason()] << std::endl;
}



template <typename T>
SNESConvergedReason PetscNonlinearSolver<T>::get_converged_reason()
{
  int ierr=0;
  
  if (this->initialized())
    {
      ierr = SNESGetConvergedReason(_snes, &_reason);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  return _reason;
}
  



//------------------------------------------------------------------
// Explicit instantiations
template class PetscNonlinearSolver<Number>;

} // namespace libMesh
 


#endif // #ifdef LIBMESH_HAVE_PETSC
