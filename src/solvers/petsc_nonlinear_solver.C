// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_mffd_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/preconditioner.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petscdmlibmesh.h"
#include "libmesh/petsc_mffd_matrix.h"

#ifdef LIBMESH_HAVE_PETSC_HYPRE
#include <HYPRE_utilities.h>
#endif

namespace libMesh
{
class ResidualContext
{
public:
  ResidualContext(PetscNonlinearSolver<Number> * solver_in, NonlinearImplicitSystem & sys_in) :
      solver(solver_in),
      sys(sys_in)
    {}

  PetscNonlinearSolver<Number> * solver;
  NonlinearImplicitSystem & sys;
};

ResidualContext
libmesh_petsc_snes_residual_helper (SNES snes, Vec x, void * ctx)
{
  LOG_SCOPE("residual()", "PetscNonlinearSolver");

  libmesh_assert(x);
  libmesh_assert(ctx);

  // No way to safety-check this cast, since we got a void *...
  PetscNonlinearSolver<Number> * solver =
    static_cast<PetscNonlinearSolver<Number> *> (ctx);

  libmesh_parallel_only(solver->comm());

  // Get the current iteration number from the snes object,
  // store it in the PetscNonlinearSolver object for possible use
  // by the user's residual function.
  {
    PetscInt n_iterations = 0;
    LibmeshPetscCall2(solver->comm(), SNESGetIterationNumber(snes, &n_iterations));
    solver->_current_nonlinear_iteration_number = cast_int<unsigned>(n_iterations);
  }

  NonlinearImplicitSystem & sys = solver->system();

  PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());

  PetscVector<Number> X_global(x, sys.comm());

  // Use the system's update() to get a good local version of the
  // parallel solution.  This operation does not modify the incoming
  // "x" vector, it only localizes information from "x" into
  // sys.current_local_solution.
  X_global.swap(X_sys);
  sys.update();
  X_global.swap(X_sys);

  // Enforce constraints (if any) exactly on the
  // current_local_solution.  This is the solution vector that is
  // actually used in the computation of the residual below, and is
  // not locked by debug-enabled PETSc the way that "x" is.
  if (solver->_exact_constraint_enforcement)
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

  return ResidualContext(solver, sys);
}

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // -----------------------------------------------------------------
  // this function monitors the nonlinear solve and checks to see
  // if we want to recalculate the preconditioner.  It only gets
  // added to the SNES instance if we're reusing the preconditioner
  PetscErrorCode
  libmesh_petsc_recalculate_monitor(SNES snes, PetscInt, PetscReal, void* ctx)
  {
    PetscFunctionBegin;

    // No way to safety-check this cast, since we got a void *...
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (ctx);

    KSP ksp;
    LibmeshPetscCall2(solver->comm(), SNESGetKSP(snes, &ksp));

    PetscInt niter;
    LibmeshPetscCall2(solver->comm(), KSPGetIterationNumber(ksp, &niter));

    if (niter > cast_int<PetscInt>(solver->reuse_preconditioner_max_linear_its()))
    {
      // -2 is a magic number for "recalculate next time you need it
      // and then not again"
      LibmeshPetscCall2(solver->comm(), SNESSetLagPreconditioner(snes, -2));
    }
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

  //-------------------------------------------------------------------
  // this function is called by PETSc at the end of each nonlinear step
  PetscErrorCode
  libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *)
  {
    PetscFunctionBegin;
    libMesh::out << "  NL step "
                 << std::setw(2) << its
                 << std::scientific
                 << ", |residual|_2 = " << fnorm
                 << std::endl;
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode
  __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_monitor(nullptr, its, fnorm, nullptr));
  }
#endif


  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  libmesh_petsc_snes_residual (SNES snes, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;

    ResidualContext rc = libmesh_petsc_snes_residual_helper(snes, x, ctx);

    libmesh_parallel_only(rc.sys.comm());

    libmesh_assert(r);
    PetscVector<Number> R(r, rc.sys.comm());

    if (rc.solver->_zero_out_residual)
      R.zero();

    //-----------------------------------------------------------------------------
    // if the user has provided both function pointers and objects only the pointer
    // will be used, so catch that as an error
    libmesh_error_msg_if(rc.solver->residual && rc.solver->residual_object,
                         "ERROR: cannot specify both a function and object to compute the Residual!");

    libmesh_error_msg_if(rc.solver->matvec && rc.solver->residual_and_jacobian_object,
                         "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!");

    if (rc.solver->residual != nullptr)
      rc.solver->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else if (rc.solver->residual_object != nullptr)
      rc.solver->residual_object->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else if (rc.solver->matvec != nullptr)
      rc.solver->matvec (*rc.sys.current_local_solution.get(), &R, nullptr, rc.sys);

    else if (rc.solver->residual_and_jacobian_object != nullptr)
    {
      auto & jac = rc.sys.get_system_matrix();

      if (rc.solver->_zero_out_jacobian)
        jac.zero();

      rc.solver->residual_and_jacobian_object->residual_and_jacobian(
          *rc.sys.current_local_solution.get(), &R, &jac, rc.sys);

      jac.close();
      if (rc.solver->_exact_constraint_enforcement)
        {
          rc.sys.get_dof_map().enforce_constraints_on_jacobian(rc.sys, &jac);
          jac.close();
        }
    }

    else
      libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");


    // Synchronize PETSc x to local solution since the local solution may be changed due to the constraints
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(rc.sys.solution.get());
    PetscVector<Number> X_global(x, rc.sys.comm());

    X_global.swap(X_sys);
    rc.sys.update();
    X_global.swap(X_sys);

    R.close();

    if (rc.solver->_exact_constraint_enforcement)
      {
        rc.sys.get_dof_map().enforce_constraints_on_residual(rc.sys, &R, rc.sys.current_local_solution.get());
        R.close();
      }

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode
  __libmesh_petsc_snes_residual (SNES snes, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_residual(snes, x, r, ctx));
  }
#endif

  //-----------------------------------------------------------------------------------------
  // this function is called by PETSc to approximate the Jacobian at X via finite differences
  PetscErrorCode
  libmesh_petsc_snes_fd_residual (SNES snes, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;

    ResidualContext rc = libmesh_petsc_snes_residual_helper(snes, x, ctx);

    libmesh_parallel_only(rc.sys.comm());

    libmesh_assert(r);
    PetscVector<Number> R(r, rc.sys.comm());

    if (rc.solver->_zero_out_residual)
      R.zero();

    if (rc.solver->fd_residual_object != nullptr)
      rc.solver->fd_residual_object->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else if (rc.solver->residual_object != nullptr)
      rc.solver->residual_object->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else
      libmesh_error_msg("Error! Unable to compute residual for forming finite difference Jacobian!");

    // Synchronize PETSc x to local solution since the local solution may be changed due to the constraints
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(rc.sys.solution.get());
    PetscVector<Number> X_global(x, rc.sys.comm());

    X_global.swap(X_sys);
    rc.sys.update();
    X_global.swap(X_sys);

    R.close();

    if (rc.solver->_exact_constraint_enforcement)
      {
        rc.sys.get_dof_map().enforce_constraints_on_residual(rc.sys, &R, rc.sys.current_local_solution.get());
        R.close();
      }

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode
  __libmesh_petsc_snes_fd_residual (SNES snes, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_fd_residual(snes, x, r, ctx));
  }
#endif

  //----------------------------------------------------------------
  // this function is called by PETSc to approximate Jacobian-vector
  // products at X via finite differences
  PetscErrorCode
  libmesh_petsc_snes_mffd_residual (SNES snes, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;

    ResidualContext rc = libmesh_petsc_snes_residual_helper(snes, x, ctx);

    libmesh_parallel_only(rc.sys.comm());

    libmesh_assert(r);
    PetscVector<Number> R(r, rc.sys.comm());

    if (rc.solver->_zero_out_residual)
      R.zero();

    if (rc.solver->mffd_residual_object != nullptr)
      rc.solver->mffd_residual_object->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else if (rc.solver->residual_object != nullptr)
      rc.solver->residual_object->residual(*rc.sys.current_local_solution.get(), R, rc.sys);

    else
      libmesh_error_msg("Error! Unable to compute residual for forming finite differenced"
                        "Jacobian-vector products!");

    // Synchronize PETSc x to local solution since the local solution may be changed due to the constraints
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(rc.sys.solution.get());
    PetscVector<Number> X_global(x, rc.sys.comm());

    X_global.swap(X_sys);
    rc.sys.update();
    X_global.swap(X_sys);

    R.close();

    if (rc.solver->_exact_constraint_enforcement)
      {
        rc.sys.get_dof_map().enforce_constraints_on_residual(rc.sys, &R, rc.sys.current_local_solution.get());
        R.close();
      }

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

  //----------------------------------------------------------
  // this function serves an interface between the petsc layer
  // and the actual mffd residual computing routine
  PetscErrorCode
  libmesh_petsc_snes_mffd_interface (void * ctx, Vec x, Vec r)
  {
    PetscFunctionBegin;

    // No way to safety-check this cast, since we got a void *...
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (ctx);

    LibmeshPetscCall2(solver->comm(), libmesh_petsc_snes_mffd_residual(solver->snes(), x, r, ctx));

#if !PETSC_VERSION_LESS_THAN(3,8,4)
#ifndef NDEBUG

    // When the user requested to reuse the nonlinear residual as the base for doing matrix-free
    // approximation of the Jacobian, we'll do a sanity check to make sure that that was safe to do
    if (solver->snes_mf_reuse_base() && (solver->comm().size() == 1) && (libMesh::n_threads() == 1))
    {
      SNES snes = solver->snes();

      KSP ksp;
      LibmeshPetscCall2(solver->comm(), SNESGetKSP(snes, &ksp));

      PetscInt ksp_it;
      LibmeshPetscCall2(solver->comm(), KSPGetIterationNumber(ksp, &ksp_it));

      SNESType snes_type;
      LibmeshPetscCall2(solver->comm(), SNESGetType(snes, &snes_type));

      libmesh_assert_msg(snes_type, "We're being called from SNES; snes_type should be non-null");

      Mat J;
      LibmeshPetscCall2(solver->comm(), SNESGetJacobian(snes, &J, NULL, NULL, NULL));
      libmesh_assert_msg(J, "We're being called from SNES; J should be non-null");

      MatType mat_type;
      LibmeshPetscCall2(solver->comm(), MatGetType(J, &mat_type));
      libmesh_assert_msg(mat_type, "We're being called from SNES; mat_type should be non-null");

      bool is_operator_mffd = strcmp(mat_type, MATMFFD) == 0;

      if ((ksp_it == PetscInt(0)) && is_operator_mffd)
      {
        bool computing_base_vector = solver->computing_base_vector();

        if (computing_base_vector)
        {
          Vec nonlinear_residual;

          LibmeshPetscCall2(solver->comm(), SNESGetFunction(snes, &nonlinear_residual, NULL, NULL));

          PetscBool vecs_equal;
          LibmeshPetscCall2(solver->comm(), VecEqual(r, nonlinear_residual, &vecs_equal));

          libmesh_error_msg_if(!(vecs_equal == PETSC_TRUE),
                               "You requested to reuse the nonlinear residual vector as the base vector for "
                               "computing the action of the matrix-free Jacobian, but the vectors are not "
                               "the same. Your physics must have states; either remove the states "
                               "from your code or make sure that you set_mf_reuse_base(false)");
        }

        // There are always exactly two function evaluations for the zeroth ksp iteration when doing
        // matrix-free approximation of the Jacobian action: one corresponding to the evaluation of
        // the base vector, and the other corresponding to evaluation of the perturbed vector. So we
        // toggle back and forth between states
        solver->set_computing_base_vector(!computing_base_vector);
      }
    }
#endif
#endif

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode
  __libmesh_petsc_snes_mffd_interface (void * ctx, Vec x, Vec r)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_mffd_interface(ctx, x, r));
  }
#endif

  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  libmesh_petsc_snes_jacobian(SNES snes, Vec x, Mat jac, Mat pc, void * ctx)
  {
    PetscFunctionBegin;

    LOG_SCOPE("jacobian()", "PetscNonlinearSolver");

    libmesh_assert(ctx);

    // No way to safety-check this cast, since we got a void *...
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (ctx);

    libmesh_parallel_only(solver->comm());

    // Get the current iteration number from the snes object,
    // store it in the PetscNonlinearSolver object for possible use
    // by the user's Jacobian function.
    {
      PetscInt n_iterations = 0;
      LibmeshPetscCall2(solver->comm(), SNESGetIterationNumber(snes, &n_iterations));
      solver->_current_nonlinear_iteration_number = cast_int<unsigned>(n_iterations);
    }

    //-----------------------------------------------------------------------------
    // if the user has provided both function pointers and objects only the pointer
    // will be used, so catch that as an error
    libmesh_error_msg_if(solver->jacobian && solver->jacobian_object,
                         "ERROR: cannot specify both a function and object to compute the Jacobian!");

    libmesh_error_msg_if(solver->matvec && solver->residual_and_jacobian_object,
                         "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!");

    NonlinearImplicitSystem & sys = solver->system();

    PetscMatrixBase<Number> * const PC = pc ? PetscMatrixBase<Number>::get_context(pc, sys.comm()) : nullptr;
    PetscMatrixBase<Number> * Jac = jac ? PetscMatrixBase<Number>::get_context(jac, sys.comm()) : nullptr;
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X_global(x, sys.comm());

    PetscMFFDMatrix<Number> mffd_jac(sys.comm());
    PetscBool p_is_shell = PETSC_FALSE;
    PetscBool j_is_mffd = PETSC_FALSE;
    PetscBool j_is_shell = PETSC_FALSE;
    if (pc)
      LibmeshPetscCall2(sys.comm(), PetscObjectTypeCompare((PetscObject)pc, MATSHELL, &p_is_shell));
    libmesh_assert(jac);
    LibmeshPetscCall2(sys.comm(), PetscObjectTypeCompare((PetscObject)jac, MATMFFD, &j_is_mffd));
    LibmeshPetscCall2(sys.comm(), PetscObjectTypeCompare((PetscObject)jac, MATSHELL, &j_is_shell));
    if (j_is_mffd == PETSC_TRUE)
      {
        libmesh_assert(!Jac);
        Jac = &mffd_jac;
        mffd_jac = jac;
      }

    // We already computed the Jacobian during the residual evaluation
    if (solver->residual_and_jacobian_object)
    {
      // We could be doing matrix-free in which case we cannot rely on closing of explicit matrices
      // that occurs during the PETSc residual callback
      if ((j_is_shell == PETSC_TRUE) || (j_is_mffd == PETSC_TRUE))
        Jac->close();

      if (pc && (p_is_shell == PETSC_TRUE))
        PC->close();

      PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
    }

    // Set the dof maps
    PC->attach_dof_map(sys.get_dof_map());
    Jac->attach_dof_map(sys.get_dof_map());

    // Use the systems update() to get a good local version of the parallel solution
    X_global.swap(X_sys);
    sys.update();
    X_global.swap(X_sys);

    // Enforce constraints (if any) exactly on the
    // current_local_solution.  This is the solution vector that is
    // actually used in the computation of the residual below, and is
    // not locked by debug-enabled PETSc the way that "x" is.
    if (solver->_exact_constraint_enforcement)
      sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    if (solver->_zero_out_jacobian)
      PC->zero();


    if (solver->jacobian != nullptr)
      solver->jacobian(*sys.current_local_solution.get(), *PC, sys);

    else if (solver->jacobian_object != nullptr)
      solver->jacobian_object->jacobian(*sys.current_local_solution.get(), *PC, sys);

    else if (solver->matvec != nullptr)
      solver->matvec(*sys.current_local_solution.get(), nullptr, PC, sys);

    else
      libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");

    PC->close();
    if (solver->_exact_constraint_enforcement)
      {
        sys.get_dof_map().enforce_constraints_on_jacobian(sys, PC);
        PC->close();
      }

    if (Jac != PC)
      {
        // Assume that shells know what they're doing
        libmesh_assert(!solver->_exact_constraint_enforcement || (j_is_mffd == PETSC_TRUE) ||
                       (j_is_shell == PETSC_TRUE));
        Jac->close();
      }

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

  // This function gets called by PETSc in place of the standard Petsc line searches
  // if a linesearch object is supplied to the PetscNonlinearSolver class. It wraps
  // the linesearch algorithm implemented on the linesearch object.
  // * "linesearch" is an object that can be used to access the non-linear and linear solution
  // vectors as well as the residual and SNES object
  // * "ctx" is the PetscNonlinearSolver context
  PetscErrorCode libmesh_petsc_linesearch_shellfunc (SNESLineSearch linesearch, void * ctx)
  {
    PetscFunctionBegin;

    // No way to safety-check this cast, since we got a void *...
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (ctx);

    libmesh_parallel_only(solver->comm());

    solver->linesearch_object->linesearch(linesearch);
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode
  __libmesh_petsc_snes_jacobian(SNES snes, Vec x, Mat jac, Mat pc, void * ctx)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_jacobian(snes, x, jac, pc, ctx));
  }
#endif

  // This function gets called by PETSc after the SNES linesearch is
  // complete.  We use it to exactly enforce any constraints on the
  // solution which may have drifted during the linear solve.  In the
  // PETSc nomenclature:
  // * "x" is the old solution vector,
  // * "y" is the search direction (Newton step) vector,
  // * "w" is the candidate solution vector, and
  // the user is responsible for setting changed_y and changed_w
  // appropriately, depending on whether or not the search
  // direction or solution vector was changed, respectively.
  PetscErrorCode libmesh_petsc_snes_postcheck(SNESLineSearch, Vec x, Vec y, Vec w, PetscBool * changed_y, PetscBool * changed_w, void * context)
  {
    PetscFunctionBegin;

    LOG_SCOPE("postcheck()", "PetscNonlinearSolver");

    // PETSc almost certainly initializes these to false already, but
    // it doesn't hurt to be explicit.
    *changed_w = PETSC_FALSE;
    *changed_y = PETSC_FALSE;

    libmesh_assert(context);

    // Cast the context to a NonlinearSolver object.
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (context);

    libmesh_parallel_only(solver->comm());

    // If the user has provided both postcheck function pointer and
    // object, this is ambiguous, so throw an error.
    libmesh_error_msg_if(solver->postcheck && solver->postcheck_object,
                         "ERROR: cannot specify both a function and object for performing the solve postcheck!");

    // It's also possible that we don't need to do anything at all, in
    // that case return early...
    NonlinearImplicitSystem & sys = solver->system();

    if (!solver->postcheck && !solver->postcheck_object)
      PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);

    // We definitely need to wrap at least "w"
    PetscVector<Number> petsc_w(w, sys.comm());

    // The user sets these flags in his/her postcheck function to
    // indicate whether they changed something.
    bool
      changed_search_direction = false,
      changed_new_soln = false;

    if (solver->postcheck || solver->postcheck_object)
      {
        PetscVector<Number> petsc_x(x, sys.comm());
        PetscVector<Number> petsc_y(y, sys.comm());

        if (solver->postcheck)
          solver->postcheck(petsc_x,
                            petsc_y,
                            petsc_w,
                            changed_search_direction,
                            changed_new_soln,
                            sys);

        else if (solver->postcheck_object)
          solver->postcheck_object->postcheck(petsc_x,
                                              petsc_y,
                                              petsc_w,
                                              changed_search_direction,
                                              changed_new_soln,
                                              sys);
      }

    // Record whether the user changed the solution or the search direction.
    if (changed_search_direction)
      *changed_y = PETSC_TRUE;

    if (changed_new_soln)
      *changed_w = PETSC_TRUE;

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  PetscErrorCode __libmesh_petsc_snes_postcheck(SNESLineSearch, Vec x, Vec y, Vec w, PetscBool * changed_y, PetscBool * changed_w, void * context)
  {
    PetscFunctionBegin;
    libmesh_deprecated();
    PetscFunctionReturn(libmesh_petsc_snes_postcheck(nullptr, x, y, w, changed_y, changed_w, context));
  }
#endif

  PetscErrorCode libmesh_petsc_snes_precheck(SNESLineSearch, Vec X, Vec Y, PetscBool * changed, void * context)
  {
    PetscFunctionBegin;

    LOG_SCOPE("precheck()", "PetscNonlinearSolver");

    // PETSc almost certainly initializes these to false already, but
    // it doesn't hurt to be explicit.
    *changed = PETSC_FALSE;

    libmesh_assert(context);

    // Cast the context to a NonlinearSolver object.
    PetscNonlinearSolver<Number> * solver =
      static_cast<PetscNonlinearSolver<Number> *> (context);

    libmesh_parallel_only(solver->comm());

    // It's possible that we don't need to do anything at all, in
    // that case return early...
    if (!solver->precheck_object)
      PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);

    // The user sets these flags in his/her postcheck function to
    // indicate whether they changed something.
    bool
      petsc_changed = false;

    auto & sys = solver->system();
    auto & x_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> petsc_x(X, sys.comm());
    PetscVector<Number> petsc_y(Y, sys.comm());

    // Use the systems update() to get a good local version of the parallel solution
    petsc_x.swap(x_sys);
    sys.update();
    petsc_x.swap(x_sys);

    // Enforce constraints (if any) exactly on the
    // current_local_solution.  This is the solution vector that is
    // actually used in the computation of residuals and Jacobians, and is
    // not locked by debug-enabled PETSc the way that "x" is.
    libmesh_assert(sys.current_local_solution.get());
    auto & local_soln = *sys.current_local_solution.get();
    if (solver->_exact_constraint_enforcement)
      sys.get_dof_map().enforce_constraints_exactly(sys, &local_soln);

    solver->precheck_object->precheck(local_soln,
                                      petsc_y,
                                      petsc_changed,
                                      sys);

    // Record whether the user changed the solution or the search direction.
    if (petsc_changed)
      *changed = PETSC_TRUE;

    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

} // end extern "C"



//---------------------------------------------------------------------
// PetscNonlinearSolver<> methods
template <typename T>
PetscNonlinearSolver<T>::PetscNonlinearSolver (sys_type & system_in) :
  NonlinearSolver<T>(system_in),
  linesearch_object(nullptr),
  _reason(SNES_CONVERGED_ITERATING/*==0*/), // Arbitrary initial value...
  _n_linear_iterations(0),
  _current_nonlinear_iteration_number(0),
  _zero_out_residual(true),
  _zero_out_jacobian(true),
  _default_monitor(true),
  _snesmf_reuse_base(true),
  _computing_base_vector(true),
  _setup_reuse(false),
  _mffd_jac(this->_communicator)
{
}



template <typename T>
PetscNonlinearSolver<T>::~PetscNonlinearSolver () = default;



template <typename T>
void PetscNonlinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      // If we don't need the preconditioner next time
      // retain the original behavior of clearing the data
      // between solves.
      if (!(reuse_preconditioner()))
        {
        // SNESReset really ought to work but replacing destroy() with
        // SNESReset causes a very slight change in behavior that
        // manifests as two failed MOOSE tests...
        _snes.destroy();
        }

      // Reset the nonlinear iteration counter.  This information is only relevant
      // *during* the solve().  After the solve is completed it should return to
      // the default value of 0.
      _current_nonlinear_iteration_number = 0;
    }
}

template <typename T>
void PetscNonlinearSolver<T>::init (const char * name)
{
  parallel_object_only();

  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      // Make only if we don't already have a retained snes
      // hanging around from the last solve
      if (!_snes)
        LibmeshPetscCall(SNESCreate(this->comm().get(), _snes.get()));

      // I believe all of the following can be safely repeated
      // even on an old snes instance from the last solve

      if (name)
        {
          libmesh_assert(std::string(name).front() != '-');
          libmesh_assert(std::string(name).back() == '_');
          LibmeshPetscCall(SNESSetOptionsPrefix(_snes, name));
        }

      // Attaching a DM to SNES.
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
      bool use_petsc_dm = libMesh::on_command_line(
          "--" + (name ? std::string(name) : std::string("")) + "use_petsc_dm");

      // This needs to be called before SNESSetFromOptions
      if (use_petsc_dm)
        this->_dm_wrapper.init_and_attach_petscdm(this->system(), _snes);
      else
#endif
      {
        WrappedPetsc<DM> dm;
        LibmeshPetscCall(DMCreate(this->comm().get(), dm.get()));
        LibmeshPetscCall(DMSetType(dm, DMLIBMESH));
        LibmeshPetscCall(DMlibMeshSetSystem(dm, this->system()));

        if (name)
          LibmeshPetscCall(DMSetOptionsPrefix(dm, name));

        LibmeshPetscCall(DMSetFromOptions(dm));
        LibmeshPetscCall(DMSetUp(dm));
        LibmeshPetscCall(SNESSetDM(_snes, dm));
        // SNES now owns the reference to dm.
      }

      setup_default_monitor();

      // If the SolverConfiguration object is provided, use it to set
      // options during solver initialization.
      if (this->_solver_configuration)
        {
          this->_solver_configuration->set_options_during_init();
        }

      if (this->_preconditioner)
        {
          KSP ksp;
          LibmeshPetscCall(SNESGetKSP (_snes, &ksp));
          PC pc;
          LibmeshPetscCall(KSPGetPC(ksp,&pc));

          this->_preconditioner->init();

          LibmeshPetscCall(PCSetType(pc, PCSHELL));
          LibmeshPetscCall(PCShellSetContext(pc,(void *)this->_preconditioner));

          //Re-Use the shell functions from petsc_linear_solver
          LibmeshPetscCall(PCShellSetSetUp(pc,libmesh_petsc_preconditioner_setup));
          LibmeshPetscCall(PCShellSetApply(pc,libmesh_petsc_preconditioner_apply));
        }
    }


  // Tell PETSc about our linesearch "post-check" function, but only
  // if the user has provided one.  There seem to be extra,
  // unnecessary residual calculations if a postcheck function is
  // attached for no reason.
  if (this->postcheck || this->postcheck_object)
    {
      SNESLineSearch linesearch;
      LibmeshPetscCall(SNESGetLineSearch(_snes, &linesearch));

      LibmeshPetscCall(SNESLineSearchSetPostCheck(linesearch, libmesh_petsc_snes_postcheck, this));
    }

  if (this->precheck_object)
    {
      SNESLineSearch linesearch;
      LibmeshPetscCall(SNESGetLineSearch(_snes, &linesearch));

      LibmeshPetscCall(SNESLineSearchSetPreCheck(linesearch, libmesh_petsc_snes_precheck, this));
    }
}


template <typename T>
SNES PetscNonlinearSolver<T>::snes(const char * name)
{
  this->init(name);
  return _snes;
}



template <typename T>
void
PetscNonlinearSolver<T>::build_mat_null_space(NonlinearImplicitSystem::ComputeVectorSubspace * computeSubspaceObject,
                                              void (*computeSubspace)(std::vector<NumericVector<Number> *> &, sys_type &),
                                              MatNullSpace * msp)
{
  parallel_object_only();

  std::vector<NumericVector<Number> *> sp;
  if (computeSubspaceObject)
    (*computeSubspaceObject)(sp, this->system());
  else
    (*computeSubspace)(sp, this->system());

  *msp = LIBMESH_PETSC_NULLPTR;
  if (sp.size())
    {
      PetscInt nmodes = cast_int<PetscInt>(sp.size());

      std::vector<Vec> modes(nmodes);
      std::vector<PetscScalar> dots(nmodes);

      for (PetscInt i=0; i<nmodes; ++i)
        {
          auto pv = cast_ptr<PetscVector<T> *>(sp[i]);

          LibmeshPetscCall(VecDuplicate(pv->vec(), &modes[i]));

          LibmeshPetscCall(VecCopy(pv->vec(), modes[i]));
        }

      // Normalize.
      LibmeshPetscCall(VecNormalize(modes[0], LIBMESH_PETSC_NULLPTR));

      for (PetscInt i=1; i<nmodes; i++)
        {
          // Orthonormalize vec[i] against vec[0:i-1]
          LibmeshPetscCall(VecMDot(modes[i], i, modes.data(), dots.data()));

          for (PetscInt j=0; j<i; j++)
            dots[j] *= -1.;

          LibmeshPetscCall(VecMAXPY(modes[i], i, dots.data(), modes.data()));

          LibmeshPetscCall(VecNormalize(modes[i], LIBMESH_PETSC_NULLPTR));
        }

      LibmeshPetscCall(MatNullSpaceCreate(this->comm().get(), PETSC_FALSE, nmodes, modes.data(), msp));

      for (PetscInt i=0; i<nmodes; ++i)
        LibmeshPetscCall(VecDestroy(&modes[i]));
    }
}

template <typename T>
std::pair<unsigned int, Real>
PetscNonlinearSolver<T>::solve (SparseMatrix<T> &  pre_in,  // System Preconditioning Matrix
                                NumericVector<T> & x_in,    // Solution vector
                                NumericVector<T> & r_in,    // Residual vector
                                const double,              // Stopping tolerance
                                const unsigned int)
{
  parallel_object_only();

  LOG_SCOPE("solve()", "PetscNonlinearSolver");
  this->init ();

  // Make sure the data passed in are really of Petsc types
  PetscMatrixBase<T> * pre = cast_ptr<PetscMatrixBase<T> *>(&pre_in);
  PetscVector<T> * x   = cast_ptr<PetscVector<T> *>(&x_in);
  PetscVector<T> * r   = cast_ptr<PetscVector<T> *>(&r_in);

  PetscInt n_iterations =0;
  // Should actually be a PetscReal, but I don't know which version of PETSc first introduced PetscReal
  Real final_residual_norm=0.;

  // We don't want to do this twice because it resets
  // SNESSetLagPreconditioner
  if ((reuse_preconditioner()) && (!_setup_reuse))
    {
      _setup_reuse = true;
      LibmeshPetscCall(SNESSetLagPreconditionerPersists(_snes, PETSC_TRUE));
      // According to the PETSC 3.16.5 docs -2 is a magic number which
      // means "recalculate the next time you need it and then not again"
      LibmeshPetscCall(SNESSetLagPreconditioner(_snes, -2));
      // Add in our callback which will trigger recalculating
      // the preconditioner when we hit reuse_preconditioner_max_linear_its
      LibmeshPetscCall(SNESMonitorSet(_snes, &libmesh_petsc_recalculate_monitor,
                                      this,
                                      NULL));
    }
  else if (!(reuse_preconditioner()))
    // This covers the case where it was enabled but was then disabled
    {
      LibmeshPetscCall(SNESSetLagPreconditionerPersists(_snes, PETSC_FALSE));
      if (_setup_reuse)
        {
          _setup_reuse = false;
          LibmeshPetscCall(SNESMonitorCancel(_snes));
          // Readd default monitor
          setup_default_monitor();
        }
    }

  LibmeshPetscCall(SNESSetFunction (_snes, r->vec(), libmesh_petsc_snes_residual, this));

  // Only set the jacobian function if we've been provided with something to call.
  // This allows a user to set their own jacobian function if they want to
  if (this->jacobian || this->jacobian_object || this->residual_and_jacobian_object)
    LibmeshPetscCall(SNESSetJacobian (_snes, pre->mat(), pre->mat(), libmesh_petsc_snes_jacobian, this));


  // Only set the nullspace if we have a way of computing it and the result is non-empty.
  if (this->nullspace || this->nullspace_object)
    {
      WrappedPetsc<MatNullSpace> msp;
      this->build_mat_null_space(this->nullspace_object, this->nullspace, msp.get());
      if (msp)
        LibmeshPetscCall(MatSetNullSpace(pre->mat(), msp));
    }

  // Only set the transpose nullspace if we have a way of computing it and the result is non-empty.
  if (this->transpose_nullspace || this->transpose_nullspace_object)
    {
#if PETSC_VERSION_LESS_THAN(3,6,0)
      libmesh_warning("MatSetTransposeNullSpace is only supported for PETSc >= 3.6, transpose nullspace will be ignored.");
#else
      WrappedPetsc<MatNullSpace> msp;
      this->build_mat_null_space(this->transpose_nullspace_object, this->transpose_nullspace, msp.get());
      if (msp)
        LibmeshPetscCall(MatSetTransposeNullSpace(pre->mat(), msp));
#endif
    }

  // Only set the nearnullspace if we have a way of computing it and the result is non-empty.
  if (this->nearnullspace || this->nearnullspace_object)
    {
      WrappedPetsc<MatNullSpace> msp;
      this->build_mat_null_space(this->nearnullspace_object, this->nearnullspace, msp.get());

      if (msp)
        LibmeshPetscCall(MatSetNearNullSpace(pre->mat(), msp));
    }

  // Have the Krylov subspace method use our good initial guess rather than 0
  KSP ksp;
  LibmeshPetscCall(SNESGetKSP (_snes, &ksp));

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values
  LibmeshPetscCall(KSPSetTolerances (ksp, this->initial_linear_tolerance, PETSC_DEFAULT,
                                     PETSC_DEFAULT, this->max_linear_iterations));

  // Set the tolerances for the non-linear solver.
  LibmeshPetscCall(SNESSetTolerances(_snes,
                                     this->absolute_residual_tolerance,
                                     this->relative_residual_tolerance,
                                     this->relative_step_tolerance,
                                     this->max_nonlinear_iterations,
                                     this->max_function_evaluations));

  // Set the divergence tolerance for the non-linear solver
#if !PETSC_VERSION_LESS_THAN(3,8,0)
  LibmeshPetscCall(SNESSetDivergenceTolerance(_snes, this->divergence_tolerance));
#endif

  //Pull in command-line options
#if PETSC_VERSION_LESS_THAN(3,7,0)
  LibmeshPetscCall(KSPSetFromOptions(ksp));
#endif
  LibmeshPetscCall(SNESSetFromOptions(_snes));

#if defined(LIBMESH_HAVE_PETSC_HYPRE) && !PETSC_VERSION_LESS_THAN(3,12,0) && defined(PETSC_HAVE_HYPRE_DEVICE)
  {
    // Make sure hypre has been initialized
    LibmeshPetscCallExternal(HYPRE_Initialize);
    PetscScalar * dummyarray;
    PetscMemType mtype;
    LibmeshPetscCall(VecGetArrayAndMemType(x->vec(), &dummyarray, &mtype));
    LibmeshPetscCall(VecRestoreArrayAndMemType(x->vec(), &dummyarray));
    if (PetscMemTypeHost(mtype))
      LibmeshPetscCallExternal(HYPRE_SetMemoryLocation, HYPRE_MEMORY_HOST);
  }
#endif

  if (this->user_presolve)
    this->user_presolve(this->system());

  //Set the preconditioning matrix
  if (this->_preconditioner)
    {
      this->_preconditioner->set_matrix(pre_in);
      this->_preconditioner->init();
    }

  // If the SolverConfiguration object is provided, use it to override
  // solver options.
  if (this->_solver_configuration)
    this->_solver_configuration->configure_solver();

  // In PETSc versions before 3.5.0, it is not possible to call
  // SNESSetUp() before the solution and rhs vectors are initialized, as
  // this triggers the
  //
  // "Solution vector cannot be right hand side vector!"
  //
  // error message. It is also not possible to call SNESSetSolution()
  // in those versions of PETSc to work around the problem, since that
  // API was removed in 3.0.0 and only restored in 3.6.0. The
  // overzealous check was moved out of SNESSetUp in PETSc 3.5.0
  // (petsc/petsc@154060b), so this code block should be safe to use
  // in 3.5.0 and later.
#if !PETSC_VERSION_LESS_THAN(3,6,0)
  LibmeshPetscCall(SNESSetSolution(_snes, x->vec()));
#endif
  LibmeshPetscCall(SNESSetUp(_snes));

  Mat J;
  LibmeshPetscCall(SNESGetJacobian(_snes, &J, LIBMESH_PETSC_NULLPTR,
                                   LIBMESH_PETSC_NULLPTR,
                                   LIBMESH_PETSC_NULLPTR));
  LibmeshPetscCall(MatMFFDSetFunction(J, libmesh_petsc_snes_mffd_interface, this));
#if !PETSC_VERSION_LESS_THAN(3,8,4)
#ifndef NDEBUG
  // If we're in debug mode, do not reuse the nonlinear function evaluation as the base for doing
  // matrix-free approximations of the Jacobian action. Instead if the user requested that we reuse
  // the base, we'll check the base function evaluation and compare it to the nonlinear residual
  // evaluation. If they are different, then we'll error and inform the user that it's unsafe to
  // reuse the base
  LibmeshPetscCall(MatSNESMFSetReuseBase(J, PETSC_FALSE));
#else
  // Resue the residual vector from SNES
  LibmeshPetscCall(MatSNESMFSetReuseBase(J, static_cast<PetscBool>(_snesmf_reuse_base)));
#endif
#endif

  SNESLineSearch linesearch;
  if (linesearch_object)
  {
    LibmeshPetscCall(SNESGetLineSearch(_snes, &linesearch));
    LibmeshPetscCall(SNESLineSearchSetType(linesearch, SNESLINESEARCHSHELL));
#if PETSC_RELEASE_GREATER_EQUALS(3, 21, 0)
    LibmeshPetscCall(SNESLineSearchShellSetApply(linesearch, libmesh_petsc_linesearch_shellfunc, this));
#else
    LibmeshPetscCall(SNESLineSearchShellSetUserFunc(linesearch, libmesh_petsc_linesearch_shellfunc, this));
#endif
  }

  LibmeshPetscCall(SNESSolve (_snes, LIBMESH_PETSC_NULLPTR, x->vec()));

  LibmeshPetscCall(SNESGetIterationNumber(_snes, &n_iterations));

  LibmeshPetscCall(SNESGetLinearSolveIterations(_snes, &_n_linear_iterations));

  // SNESGetFunction has been around forever and should work on all
  // versions of PETSc.  This is also now the recommended approach
  // according to the documentation for the PETSc 3.5.1 release:
  // http://www.mcs.anl.gov/petsc/documentation/changes/35.html
  Vec f;
  LibmeshPetscCall(SNESGetFunction(_snes, &f, 0, 0));
  LibmeshPetscCall(VecNorm(f, NORM_2, pPR(&final_residual_norm)));

  // Get and store the reason for convergence
  LibmeshPetscCall(SNESGetConvergedReason(_snes, &_reason));

  //Based on Petsc 2.3.3 documentation all diverged reasons are negative
  this->converged = (_reason >= 0);

  // Reset data structure
  this->clear();

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
  if (this->initialized())
    LibmeshPetscCall(SNESGetConvergedReason(_snes, &_reason));

  return _reason;
}

template <typename T>
int PetscNonlinearSolver<T>::get_total_linear_iterations()
{
  return _n_linear_iterations;
}

template <typename T>
void PetscNonlinearSolver<T>::setup_default_monitor()
{
  if (_default_monitor)
    LibmeshPetscCall(
      SNESMonitorSet(_snes, libmesh_petsc_snes_monitor, this, LIBMESH_PETSC_NULLPTR));
}

template <typename T>
bool PetscNonlinearSolver<T>::reuse_preconditioner() const
{
  return this->_reuse_preconditioner;
}

template <typename T>
unsigned int PetscNonlinearSolver<T>::reuse_preconditioner_max_linear_its() const
{
  return this->_reuse_preconditioner_max_linear_its;
}

template <typename T>
void PetscNonlinearSolver<T>::force_new_preconditioner()
{
  // Easiest way is just to clear everything out
  this->_is_initialized = false;
  _snes.destroy();
  _setup_reuse = false;
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscNonlinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
