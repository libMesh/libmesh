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

// Local Includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/preconditioner.h"

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
    //PetscErrorCode ierr=0;

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
  __libmesh_petsc_snes_residual (SNES snes, Vec x, Vec r, void *ctx)
  {
    START_LOG("residual()", "PetscNonlinearSolver");

    PetscErrorCode ierr=0;

    libmesh_assert(x);
    libmesh_assert(r);
    libmesh_assert(ctx);

    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);

    // Get the current iteration number from the snes object,
    // store it in the PetscNonlinearSolver object for possible use
    // by the user's residual function.
    {
      PetscInt n_iterations = 0;
      ierr = SNESGetIterationNumber(snes, &n_iterations);
      CHKERRABORT(solver->comm().get(),ierr);
      solver->_current_nonlinear_iteration_number = static_cast<unsigned>(n_iterations);
    }

    NonlinearImplicitSystem &sys = solver->system();

    PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
    PetscVector<Number>& R_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.rhs);
    PetscVector<Number> X_global(x, sys.comm()), R(r, sys.comm());

    // Use the systems update() to get a good local version of the parallel solution
    X_global.swap(X_sys);
    R.swap(R_sys);

    sys.get_dof_map().enforce_constraints_exactly(sys);

    sys.update();

    //Swap back
    X_global.swap(X_sys);
    R.swap(R_sys);

    if (solver->_zero_out_residual)
      R.zero();

    //-----------------------------------------------------------------------------
    // if the user has provided both function pointers and objects only the pointer
    // will be used, so catch that as an error
    if (solver->residual && solver->residual_object)
      {
	libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Residual!" << std::endl;
	libmesh_error();
      }

    if (solver->matvec && solver->residual_and_jacobian_object)
      {
	libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
	libmesh_error();
      }
    //-----------------------------------------------------------------------------

    if      (solver->residual != NULL)                     solver->residual                                            (*sys.current_local_solution.get(), R, sys);
    else if (solver->residual_object != NULL)              solver->residual_object->residual                           (*sys.current_local_solution.get(), R, sys);
    else if (solver->matvec   != NULL)                     solver->matvec                                              (*sys.current_local_solution.get(), &R, NULL, sys);
    else if (solver->residual_and_jacobian_object != NULL) solver->residual_and_jacobian_object->residual_and_jacobian (*sys.current_local_solution.get(), &R, NULL, sys);
    else libmesh_error();

    R.close();
    X_global.close();

    STOP_LOG("residual()", "PetscNonlinearSolver");

    return ierr;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __libmesh_petsc_snes_jacobian (SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    START_LOG("jacobian()", "PetscNonlinearSolver");

    PetscErrorCode ierr=0;

    libmesh_assert(ctx);

    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);

    // Get the current iteration number from the snes object,
    // store it in the PetscNonlinearSolver object for possible use
    // by the user's Jacobian function.
    {
      PetscInt n_iterations = 0;
      ierr = SNESGetIterationNumber(snes, &n_iterations);
      CHKERRABORT(solver->comm().get(),ierr);
      solver->_current_nonlinear_iteration_number = static_cast<unsigned>(n_iterations);
    }

    NonlinearImplicitSystem &sys = solver->system();

    PetscMatrix<Number> PC(*pc, sys.comm());
    PetscMatrix<Number> Jac(*jac, sys.comm());
    PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
    PetscMatrix<Number>& Jac_sys = *libmesh_cast_ptr<PetscMatrix<Number>*>(sys.matrix);
    PetscVector<Number> X_global(x, sys.comm());

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

    if (solver->_zero_out_jacobian)
      PC.zero();

    //-----------------------------------------------------------------------------
    // if the user has provided both function pointers and objects only the pointer
    // will be used, so catch that as an error
    if (solver->jacobian && solver->jacobian_object)
      {
	libMesh::err << "ERROR: cannot specify both a function and object to compute the Jacobian!" << std::endl;
	libmesh_error();
      }

    if (solver->matvec && solver->residual_and_jacobian_object)
      {
	libMesh::err << "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!" << std::endl;
	libmesh_error();
      }
    //-----------------------------------------------------------------------------

    if      (solver->jacobian != NULL)                     solver->jacobian                                            (*sys.current_local_solution.get(), PC, sys);
    else if (solver->jacobian_object != NULL)              solver->jacobian_object->jacobian                           (*sys.current_local_solution.get(), PC, sys);
    else if (solver->matvec   != NULL)                     solver->matvec                                              (*sys.current_local_solution.get(), NULL, &PC, sys);
    else if (solver->residual_and_jacobian_object != NULL) solver->residual_and_jacobian_object->residual_and_jacobian (*sys.current_local_solution.get(), NULL, &PC, sys);
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
PetscNonlinearSolver<T>::PetscNonlinearSolver (sys_type& system_in) :
    NonlinearSolver<T>(system_in),
    _reason(SNES_CONVERGED_ITERATING/*==0*/), // Arbitrary initial value...
    _n_linear_iterations(0),
    _current_nonlinear_iteration_number(0),
    _zero_out_residual(true),
    _zero_out_jacobian(true),
    _default_monitor(true)
{
}



template <typename T>
PetscNonlinearSolver<T>::~PetscNonlinearSolver ()
{
  this->clear ();
}



template <typename T>
void PetscNonlinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      PetscErrorCode ierr=0;

      ierr = LibMeshSNESDestroy(&_snes);
             LIBMESH_CHKERRABORT(ierr);

      // Reset the nonlinear iteration counter.  This information is only relevant
      // *during* the solve().  After the solve is completed it should return to
      // the default value of 0.
      _current_nonlinear_iteration_number = 0;
    }
}



template <typename T>
void PetscNonlinearSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

#if PETSC_VERSION_LESS_THAN(2,1,2)
      // At least until Petsc 2.1.1, the SNESCreate had a different calling syntax.
      // The second argument was of type SNESProblemType, and could have a value of
      // either SNES_NONLINEAR_EQUATIONS or SNES_UNCONSTRAINED_MINIMIZATION.
      ierr = SNESCreate(this->comm().get(), SNES_NONLINEAR_EQUATIONS, &_snes);
             LIBMESH_CHKERRABORT(ierr);

#else

      ierr = SNESCreate(this->comm().get(),&_snes);
             LIBMESH_CHKERRABORT(ierr);

#endif

      if (_default_monitor)
      {
#if PETSC_VERSION_LESS_THAN(2,3,3)
        ierr = SNESSetMonitor (_snes, __libmesh_petsc_snes_monitor,
			       this, PETSC_NULL);
#else
        // API name change in PETSc 2.3.3
        ierr = SNESMonitorSet (_snes, __libmesh_petsc_snes_monitor,
			       this, PETSC_NULL);
#endif
        LIBMESH_CHKERRABORT(ierr);
      }

#if PETSC_VERSION_LESS_THAN(3,1,0)
      // Cannot call SNESSetOptions before SNESSetFunction when using
      // any matrix free options with PETSc 3.1.0+
      ierr = SNESSetFromOptions(_snes);
             LIBMESH_CHKERRABORT(ierr);
#endif

      if(this->_preconditioner)
      {
        KSP ksp;
        ierr = SNESGetKSP (_snes, &ksp);
               LIBMESH_CHKERRABORT(ierr);
        PC pc;
        ierr = KSPGetPC(ksp,&pc);
               LIBMESH_CHKERRABORT(ierr);

        this->_preconditioner->init();

        PCSetType(pc, PCSHELL);
        PCShellSetContext(pc,(void*)this->_preconditioner);

        //Re-Use the shell functions from petsc_linear_solver
        PCShellSetSetUp(pc,__libmesh_petsc_preconditioner_setup);
        PCShellSetApply(pc,__libmesh_petsc_preconditioner_apply);
      }
    }
}

#if !PETSC_VERSION_LESS_THAN(3,3,0)
template <typename T>
void
PetscNonlinearSolver<T>::build_mat_null_space(NonlinearImplicitSystem::ComputeVectorSubspace* computeSubspaceObject,
                                              void (*computeSubspace)(std::vector<NumericVector<Number>*>&, sys_type&),
                                              MatNullSpace *msp)
{
  PetscErrorCode ierr;
  std::vector<NumericVector<Number>* > sp;
  if (computeSubspaceObject)
    (*computeSubspaceObject)(sp, this->system());
  else
    (*computeSubspace)(sp, this->system());

  *msp = PETSC_NULL;
  if (sp.size())
    {
      Vec *modes;
      PetscScalar *dots;
      PetscInt nmodes = sp.size();

#if PETSC_RELEASE_LESS_THAN(3,5,0)
      ierr = PetscMalloc2(nmodes,Vec,&modes,nmodes,PetscScalar,&dots);
#else
      ierr = PetscMalloc2(nmodes,&modes,nmodes,&dots);
#endif
      LIBMESH_CHKERRABORT(ierr);

      for (PetscInt i=0; i<nmodes; ++i)
        {
          PetscVector<T>* pv = libmesh_cast_ptr<PetscVector<T>*>(sp[i]);
          Vec v = pv->vec();

          ierr = VecDuplicate(v, modes+i);
          LIBMESH_CHKERRABORT(ierr);

          ierr = VecCopy(v,modes[i]);
          LIBMESH_CHKERRABORT(ierr);
        }

      // Normalize.
      ierr = VecNormalize(modes[0],PETSC_NULL);
      LIBMESH_CHKERRABORT(ierr);

      for (PetscInt i=1; i<nmodes; i++)
        {
          // Orthonormalize vec[i] against vec[0:i-1]
          ierr = VecMDot(modes[i],i,modes,dots);
          LIBMESH_CHKERRABORT(ierr);

          for (PetscInt j=0; j<i; j++)
            dots[j] *= -1.;

          ierr = VecMAXPY(modes[i],i,dots,modes);
          LIBMESH_CHKERRABORT(ierr);

          ierr = VecNormalize(modes[i],PETSC_NULL);
          LIBMESH_CHKERRABORT(ierr);
        }

      ierr = MatNullSpaceCreate(this->comm().get(), PETSC_FALSE, nmodes, modes, msp);
      LIBMESH_CHKERRABORT(ierr);

      for (PetscInt i=0; i<nmodes; ++i)
        {
          ierr = VecDestroy(modes+i);
          LIBMESH_CHKERRABORT(ierr);
        }

      ierr = PetscFree2(modes,dots);
      LIBMESH_CHKERRABORT(ierr);
    }
}
#endif

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

  PetscErrorCode ierr=0;
  PetscInt n_iterations =0;
  // Should actually be a PetscReal, but I don't know which version of PETSc first introduced PetscReal
  Real final_residual_norm=0.;

  ierr = SNESSetFunction (_snes, r->vec(), __libmesh_petsc_snes_residual, this);
         LIBMESH_CHKERRABORT(ierr);

   // Only set the jacobian function if we've been provided with something to call.
   // This allows a user to set their own jacobian function if they want to
   if (this->jacobian || this->jacobian_object || this->residual_and_jacobian_object)
   {
     ierr = SNESSetJacobian (_snes, jac->mat(), jac->mat(), __libmesh_petsc_snes_jacobian, this);
     LIBMESH_CHKERRABORT(ierr);
   }
#if !PETSC_VERSION_LESS_THAN(3,3,0)
   // Only set the nullspace if we have a way of computing it and the result is non-empty.
   if (this->nullspace || this->nullspace_object)
   {
     MatNullSpace msp;
     this->build_mat_null_space(this->nullspace_object, this->nullspace, &msp);
     if (msp)
       {
         ierr = MatSetNullSpace(jac->mat(), msp);
         LIBMESH_CHKERRABORT(ierr);

         ierr = MatNullSpaceDestroy(&msp);
         LIBMESH_CHKERRABORT(ierr);
       }
   }

   // Only set the nearnullspace if we have a way of computing it and the result is non-empty.
   if (this->nearnullspace || this->nearnullspace_object)
   {
     MatNullSpace msp = PETSC_NULL;
     this->build_mat_null_space(this->nearnullspace_object, this->nearnullspace, &msp);

     if(msp) {
       ierr = MatSetNearNullSpace(jac->mat(), msp);
       LIBMESH_CHKERRABORT(ierr);

       ierr = MatNullSpaceDestroy(&msp);
       LIBMESH_CHKERRABORT(ierr);
     }
   }
#endif
   // Have the Krylov subspace method use our good initial guess rather than 0
   KSP ksp;
   ierr = SNESGetKSP (_snes, &ksp);
          LIBMESH_CHKERRABORT(ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values
  ierr = KSPSetTolerances (ksp, this->initial_linear_tolerance, PETSC_DEFAULT,
                           PETSC_DEFAULT, this->max_linear_iterations);
         LIBMESH_CHKERRABORT(ierr);

  // Set the tolerances for the non-linear solver.
  ierr = SNESSetTolerances(_snes, this->absolute_residual_tolerance, this->relative_residual_tolerance,
                           this->relative_step_tolerance, this->max_nonlinear_iterations, this->max_function_evaluations);
         LIBMESH_CHKERRABORT(ierr);

  //Pull in command-line options
  KSPSetFromOptions(ksp);
  SNESSetFromOptions(_snes);

  if (this->user_presolve)
    this->user_presolve(this->system());

  //Set the preconditioning matrix
  if(this->_preconditioner)
  {
    this->_preconditioner->set_matrix(jac_in);
    this->_preconditioner->init();
  }

//    ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
//           LIBMESH_CHKERRABORT(ierr);

// Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
// the last one being a pointer to an int to hold the number of iterations required.
# if PETSC_VERSION_LESS_THAN(2,2,0)

 ierr = SNESSolve (_snes, x->vec(), &n_iterations);
        LIBMESH_CHKERRABORT(ierr);

// 2.2.x style
#elif PETSC_VERSION_LESS_THAN(2,3,0)

 ierr = SNESSolve (_snes, x->vec());
        LIBMESH_CHKERRABORT(ierr);

// 2.3.x & newer style
#else

  ierr = SNESSolve (_snes, PETSC_NULL, x->vec());
         LIBMESH_CHKERRABORT(ierr);

  ierr = SNESGetIterationNumber(_snes,&n_iterations);
         LIBMESH_CHKERRABORT(ierr);

  ierr = SNESGetLinearSolveIterations(_snes, &_n_linear_iterations);
         LIBMESH_CHKERRABORT(ierr);

  ierr = SNESGetFunctionNorm(_snes,&final_residual_norm);
	 LIBMESH_CHKERRABORT(ierr);

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
  PetscErrorCode ierr=0;

  if (this->initialized())
    {
      ierr = SNESGetConvergedReason(_snes, &_reason);
      LIBMESH_CHKERRABORT(ierr);
    }

  return _reason;
}

template <typename T>
int PetscNonlinearSolver<T>::get_total_linear_iterations()
{
  return _n_linear_iterations;
}


//------------------------------------------------------------------
// Explicit instantiations
template class PetscNonlinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
