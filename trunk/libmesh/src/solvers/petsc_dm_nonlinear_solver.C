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


#include "petsc_macro.h"

/*
  This only works with a recent petsc-dev (post petsc-3.2).
  Replace with a test for petsc-3.3 after it's released.
*/
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,2,0) && !PETSC_VERSION_RELEASE


// C++ includes

// Local Includes
#include "libmesh_common.h"
#include "nonlinear_implicit_system.h"
#include "petsc_dm_nonlinear_solver.h"
#include "petsc_linear_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "dof_map.h"
#include "preconditioner.h"


EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMFunction_libMesh"
PetscErrorCode DMFunction_libMesh(DM dm, Vec x, Vec r)
{
  PetscFunctionBegin;
  libmesh_assert (x != NULL);
  libmesh_assert (r != NULL);

  NonlinearImplicitSystem* _sys;
  PetscDMGetSystem(dm, _sys);
  NonlinearImplicitSystem& sys = *_sys;
  PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>* >(sys.solution.get());
  PetscVector<Number>& R_sys = *libmesh_cast_ptr<PetscVector<Number>* >(sys.rhs);
  PetscVector<Number> X_global(x), R(r);

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);
  R.swap(R_sys);

  _sys->get_dof_map().enforce_constraints_exactly(*_sys);
  _sys->update();

  // Swap back
  X_global.swap(X_sys);
  R.swap(R_sys);
  R.zero();

  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  if (_sys->nonlinear_solver->residual && _sys->nonlinear_solver->residual_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Residual!" << std::endl;
      libmesh_error();
    }

  if (_sys->nonlinear_solver->matvec && _sys->nonlinear_solver->residual_and_jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
      libmesh_error();
    }

  if (_sys->nonlinear_solver->residual != NULL)
    _sys->nonlinear_solver->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->residual_object != NULL)
    _sys->nonlinear_solver->residual_object->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->matvec   != NULL)
    _sys->nonlinear_solver->matvec(*(_sys->current_local_solution.get()), &R, NULL, *_sys);

  else if (_sys->nonlinear_solver->residual_and_jacobian_object != NULL)
    _sys->nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(_sys->current_local_solution.get()), &R, NULL, *_sys);

  else
    libmesh_error();

  R.close();
  X_global.close();

  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMJacobian_libMesh"
PetscErrorCode DMJacobian_libMesh(DM dm, Vec x, Mat jac, Mat pc, MatStructure *msflag)
{
  PetscFunctionBegin;

  NonlinearImplicitSystem* _sys;
  PetscDMGetSystem(dm, _sys);
  NonlinearImplicitSystem& sys = *_sys;

  PetscMatrix<Number>  PC(pc);
  PetscMatrix<Number>  Jac(jac);
  PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
  PetscMatrix<Number>& Jac_sys = *libmesh_cast_ptr<PetscMatrix<Number>*>(sys.matrix);
  PetscVector<Number>  X_global(x);

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

  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  if (sys.nonlinear_solver->jacobian && sys.nonlinear_solver->jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Jacobian!" << std::endl;
      libmesh_error();
    }

  if (sys.nonlinear_solver->matvec && sys.nonlinear_solver->residual_and_jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
      libmesh_error();
    }

  if (sys.nonlinear_solver->jacobian != NULL)
    sys.nonlinear_solver->jacobian(*(sys.current_local_solution.get()), PC, sys);

  else if (sys.nonlinear_solver->jacobian_object != NULL)
    sys.nonlinear_solver->jacobian_object->jacobian(*(sys.current_local_solution.get()), PC, sys);

  else if (sys.nonlinear_solver->matvec != NULL)
    sys.nonlinear_solver->matvec(*(sys.current_local_solution.get()), NULL, &PC, sys);

  else if (sys.nonlinear_solver->residual_and_jacobian_object != NULL)
    sys.nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(sys.current_local_solution.get()), NULL, &PC, sys);

  else
    libmesh_error();

  PC.close();
  Jac.close();
  X_global.close();

  *msflag = SAME_NONZERO_PATTERN;

  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMVariableBounds_libMesh"
PetscErrorCode DMVariableBounds_libMesh(DM dm, Vec xl, Vec xu)
{
  PetscFunctionBegin;
  NonlinearImplicitSystem* _sys;
  PetscDMGetSystem(dm, _sys);
  NonlinearImplicitSystem& sys = *_sys;
  PetscVector<Number> XL(xl);
  PetscVector<Number> XU(xu);

  if (sys.nonlinear_solver->bounds != NULL)
    sys.nonlinear_solver->bounds(XL,XU,sys);
  else if (sys.nonlinear_solver->bounds_object != NULL)
    sys.nonlinear_solver->bounds_object->bounds(XL,XU, sys);
  else
    libmesh_error();

  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_libMesh"
PetscErrorCode DMCreateGlobalVector_libMesh(DM dm, Vec *x)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  NumericVector<Number>* nv = (dlm->sys->solution).get();
  PetscVector<Number>*   pv = dynamic_cast<PetscVector<Number>*>(nv);
  Vec                    v  = pv->vec();
  ierr = VecDuplicate(v,x); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMCreateMatrix_libMesh"
PetscErrorCode DMCreateMatrix_libMesh(DM dm, const MatType, Mat *A)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  *A = (dynamic_cast<PetscMatrix<Number>*>(dlm->sys->matrix))->mat();
  ierr = PetscObjectReference((PetscObject)(*A)); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMSetUp_libMesh"
PetscErrorCode  DMSetUp_libMesh(DM dm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  ierr = DMSetFunction(dm, DMFunction_libMesh); CHKERRQ(ierr);
  ierr = DMSetJacobian(dm, DMJacobian_libMesh); CHKERRQ(ierr);
  if (dlm->sys->nonlinear_solver->bounds || dlm->sys->nonlinear_solver->bounds_object)
    ierr = DMSetVariableBounds(dm, DMVariableBounds_libMesh); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMDestroy_libMesh"
PetscErrorCode  DMDestroy_libMesh(DM)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMCreate_libMesh"
PetscErrorCode  DMCreate_libMesh(DM dm)
{
  PetscErrorCode ierr;
  DM_libMesh     *dlm;

  PetscFunctionBegin;
  PetscValidPointer(dm,1);
  ierr = PetscNewLog(dm,DM_libMesh,&dlm);CHKERRQ(ierr);
  dm->data = dlm;

  dm->ops->createglobalvector = DMCreateGlobalVector_libMesh;
  dm->ops->createlocalvector  = 0; // DMCreateLocalVector_libMesh;
  dm->ops->getcoloring        = 0; // DMGetColoring_libMesh;
  dm->ops->creatematrix       = DMCreateMatrix_libMesh;
  dm->ops->createinterpolation= 0; // DMCreateInterpolation_libMesh;

  dm->ops->refine             = 0; // DMRefine_libMesh;
  dm->ops->coarsen            = 0; // DMCoarsen_libMesh;
  dm->ops->refinehierarchy    = 0; // DMRefineHierarchy_libMesh;
  dm->ops->coarsenhierarchy   = 0; // DMCoarsenHierarchy_libMesh;
  dm->ops->getinjection       = 0; // DMGetInjection_libMesh;
  dm->ops->getaggregates      = 0; // DMGetAggregates_libMesh;
  dm->ops->destroy            = DMDestroy_libMesh;
  dm->ops->view               = 0;
  dm->ops->setfromoptions     = 0; // DMSetFromOptions_libMesh;
  dm->ops->setup              = DMSetUp_libMesh;

  PetscFunctionReturn(0);
}

EXTERN_C_END





namespace libMesh {
  PetscBool PetscDMRegistered = PETSC_FALSE;
  void PetscDMRegister()
  {
    if (PetscDMRegistered)
      return;

    PetscErrorCode ierr;
    ierr = DMRegister(DMLIBMESH, PETSC_NULL, "DMCreate_libMesh", DMCreate_libMesh); CHKERRABORT(libMesh::COMM_WORLD, ierr);
    PetscDMRegistered = PETSC_TRUE;
  }

  void PetscDMSetSystem(DM dm, NonlinearImplicitSystem& sys)
  {
    DM_libMesh *dlm = (DM_libMesh *)(dm->data);
    dlm->sys =&sys;
  }

  void PetscDMGetSystem(DM dm, NonlinearImplicitSystem*& sys)
  {
    DM_libMesh *dlm = (DM_libMesh *)(dm->data);
    sys = dlm->sys;
  }

  template <typename T>
  PetscDMNonlinearSolver<T>::PetscDMNonlinearSolver(sys_type& system) :
    PetscNonlinearSolver<T>(system)
  {
    PetscDMRegister();
  }



  template <typename T>
  inline
  PetscDMNonlinearSolver<T>::~PetscDMNonlinearSolver ()
  {
    this->clear ();
  }



  template <typename T>
  void PetscDMNonlinearSolver<T>::init()
  {
    PetscErrorCode ierr;
    DM dm;
    this->PetscNonlinearSolver<T>::init();

    // Attaching a DM with the function and Jacobian callbacks to SNES.
    ierr = DMCreate(libMesh::COMM_WORLD, &dm); CHKERRABORT(libMesh::COMM_WORLD, ierr);
    ierr = DMSetType(dm, DMLIBMESH);           CHKERRABORT(libMesh::COMM_WORLD, ierr);
    PetscDMSetSystem(dm, this->system());
    ierr = DMSetUp(dm);                        CHKERRABORT(libMesh::COMM_WORLD, ierr);
    ierr = SNESSetDM(this->_snes, dm);         CHKERRABORT(libMesh::COMM_WORLD, ierr);

    // SNES now owns the reference to dm.
    ierr = DMDestroy(&dm);                     CHKERRABORT(libMesh::COMM_WORLD, ierr);

    if (this->bounds || this->bounds_object)
      {
	ierr = SNESSetType(this->_snes, SNESVIRS);
	CHKERRABORT(libMesh::COMM_WORLD, ierr);
      }

    KSP ksp;
    ierr = SNESGetKSP (this->_snes, &ksp);     CHKERRABORT(libMesh::COMM_WORLD,ierr);

    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values
    ierr = KSPSetTolerances (ksp, this->initial_linear_tolerance, PETSC_DEFAULT,PETSC_DEFAULT, this->max_linear_iterations); CHKERRABORT(libMesh::COMM_WORLD,ierr);

    // Set the tolerances for the non-linear solver.
    ierr = SNESSetTolerances(this->_snes,
			     this->absolute_residual_tolerance,
			     this->relative_residual_tolerance,
			     this->absolute_step_tolerance,
			     this->max_nonlinear_iterations,
			     this->max_function_evaluations);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    //Pull in command-line options
    KSPSetFromOptions(ksp);
    SNESSetFromOptions(this->_snes);
  }


  template <typename T>
  std::pair<unsigned int, Real>
  PetscDMNonlinearSolver<T>::solve (SparseMatrix<T>& jac_in,  // System Jacobian Matrix
				    NumericVector<T>& x_in,   // Solution vector
				    NumericVector<T>& r_in,   // Residual vector
				    const double,             // Stopping tolerance
				    const unsigned int)
  {
    START_LOG("solve()", "PetscNonlinearSolver");
    this->init ();

    // Make sure the data passed in are really of Petsc types
    libmesh_cast_ptr<PetscMatrix<T>*>(&jac_in);
    libmesh_cast_ptr<PetscVector<T>*>(&r_in);

    // Extract solution vector
    PetscVector<T>* x = libmesh_cast_ptr<PetscVector<T>*>(&x_in);

    int ierr=0;
    int n_iterations =0;

    // Should actually be a PetscReal, but I don't know which version of PETSc first introduced PetscReal
    Real final_residual_norm=0.;

    if (this->user_presolve)
      this->user_presolve(this->system());

    //Set the preconditioning matrix
    if (this->_preconditioner)
      this->_preconditioner->set_matrix(jac_in);

    ierr = SNESSolve (this->_snes, PETSC_NULL, x->vec());
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = SNESGetIterationNumber(this->_snes,&n_iterations);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = SNESGetLinearSolveIterations(this->_snes, &this->_n_linear_iterations);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = SNESGetFunctionNorm(this->_snes,&final_residual_norm);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    // Get and store the reason for convergence
    SNESGetConvergedReason(this->_snes, &this->_reason);

    //Based on Petsc 2.3.3 documentation all diverged reasons are negative
    this->converged = (this->_reason >= 0);

    this->clear();

    STOP_LOG("solve()", "PetscNonlinearSolver");

    // return the # of its. and the final residual norm.
    return std::make_pair(n_iterations, final_residual_norm);
  }


  //------------------------------------------------------------------
  // Explicit instantiations
  template class PetscDMNonlinearSolver<Number>;

} // namespace libMesh


#endif // #if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,2,0)
