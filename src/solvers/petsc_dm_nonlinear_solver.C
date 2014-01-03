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


#include "libmesh/petsc_macro.h"

// This only works with petsc-3.3 and above.
#if !PETSC_VERSION_LESS_THAN(3,3,0)

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_dm_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/preconditioner.h"

extern "C"
{
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif
}

// PETSc extern definition
EXTERN_C_BEGIN
PetscErrorCode DMCreate_libMesh(DM);
EXTERN_C_END

#include <libmesh/petscdmlibmesh.h>

namespace libMesh {
  PetscBool PetscDMRegistered = PETSC_FALSE;
  void PetscDMRegister()
  {
    if (PetscDMRegistered)
      return;

    PetscErrorCode ierr;
#if PETSC_RELEASE_LESS_THAN(3,4,0)
    ierr = DMRegister(DMLIBMESH, PETSC_NULL, "DMCreate_libMesh", DMCreate_libMesh); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
    ierr = DMRegister(DMLIBMESH, DMCreate_libMesh); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
    PetscDMRegistered = PETSC_TRUE;
  }


  template <typename T>
  PetscDMNonlinearSolver<T>::PetscDMNonlinearSolver(sys_type& system_in) :
    PetscNonlinearSolver<T>(system_in)
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
    ierr = DMCreateLibMesh(this->comm().get(), this->system(), &dm); LIBMESH_CHKERRABORT(ierr);
    ierr = DMSetFromOptions(dm);               LIBMESH_CHKERRABORT(ierr);
    ierr = DMSetUp(dm);                        LIBMESH_CHKERRABORT(ierr);
    ierr = SNESSetDM(this->_snes, dm);         LIBMESH_CHKERRABORT(ierr);
    // SNES now owns the reference to dm.
    ierr = DMDestroy(&dm);                     LIBMESH_CHKERRABORT(ierr);
    KSP ksp;
    ierr = SNESGetKSP (this->_snes, &ksp);     LIBMESH_CHKERRABORT(ierr);

    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values
    ierr = KSPSetTolerances (ksp, this->initial_linear_tolerance, PETSC_DEFAULT,PETSC_DEFAULT, this->max_linear_iterations); LIBMESH_CHKERRABORT(ierr);

    // Set the tolerances for the non-linear solver.
    ierr = SNESSetTolerances(this->_snes,
			     this->absolute_residual_tolerance,
			     this->relative_residual_tolerance,
			     this->absolute_step_tolerance,
			     this->max_nonlinear_iterations,
			     this->max_function_evaluations);
    LIBMESH_CHKERRABORT(ierr);

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

    PetscErrorCode ierr=0;
    PetscInt n_iterations =0;


    PetscReal final_residual_norm=0.;

    if (this->user_presolve)
      this->user_presolve(this->system());

    //Set the preconditioning matrix
    if (this->_preconditioner)
      this->_preconditioner->set_matrix(jac_in);

    ierr = SNESSolve (this->_snes, PETSC_NULL, x->vec());
    LIBMESH_CHKERRABORT(ierr);

    ierr = SNESGetIterationNumber(this->_snes,&n_iterations);
    LIBMESH_CHKERRABORT(ierr);

    ierr = SNESGetLinearSolveIterations(this->_snes, &this->_n_linear_iterations);
    LIBMESH_CHKERRABORT(ierr);

    ierr = SNESGetFunctionNorm(this->_snes,&final_residual_norm);
    LIBMESH_CHKERRABORT(ierr);

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


#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
