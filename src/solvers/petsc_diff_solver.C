// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_diff_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_auto_fieldsplit.h"

#ifdef LIBMESH_HAVE_PETSC

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
  // Function to hand to PETSc's SNES,
  // which monitors convergence at X
  PetscErrorCode
  __libmesh_petsc_diff_solver_monitor (SNES snes,
                                       PetscInt its,
                                       PetscReal fnorm,
                                       void * ctx)
  {
    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver *> (ctx));

    if (solver.verbose)
      libMesh::out << "  PetscDiffSolver step " << its
                   << ", |residual|_2 = " << fnorm << std::endl;
    if (solver.linear_solution_monitor.get()) {
      int ierr = 0;

      Vec petsc_delta_u;
      ierr = SNESGetSolutionUpdate(snes, &petsc_delta_u);
      CHKERRABORT(solver.comm().get(), ierr);
      PetscVector<Number> delta_u(petsc_delta_u, solver.comm());
      delta_u.close();

      Vec petsc_u;
      ierr = SNESGetSolution(snes, &petsc_u);
      CHKERRABORT(solver.comm().get(), ierr);
      PetscVector<Number> u(petsc_u, solver.comm());
      u.close();

      Vec petsc_res;
      ierr = SNESGetFunction(snes, &petsc_res, libmesh_nullptr, libmesh_nullptr);
      CHKERRABORT(solver.comm().get(), ierr);
      PetscVector<Number> res(petsc_res, solver.comm());
      res.close();

      (*solver.linear_solution_monitor)(
                                        delta_u, delta_u.l2_norm(),
                                        u, u.l2_norm(),
                                        res, res.l2_norm(), its);
    }
    return 0;
  }

  // Functions to hand to PETSc's SNES,
  // which compute the residual or jacobian at X
  PetscErrorCode
  __libmesh_petsc_diff_solver_residual (SNES, Vec x, Vec r, void * ctx)
  {
    libmesh_assert(x);
    libmesh_assert(r);
    libmesh_assert(ctx);

    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver*> (ctx));
    ImplicitSystem & sys = solver.system();

    if (solver.verbose)
      libMesh::out << "Assembling the residual" << std::endl;

    PetscVector<Number> & X_system =
      *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> & R_system =
      *cast_ptr<PetscVector<Number> *>(sys.rhs);
    PetscVector<Number> X_input(x, sys.comm()), R_input(r, sys.comm());

    // DiffSystem assembles from the solution and into the rhs, so swap
    // those with our input vectors before assembling.  They'll probably
    // already be references to the same vectors, but PETSc might do
    // something tricky.
    X_input.swap(X_system);
    R_input.swap(R_system);

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(true, false);
    R_system.close();

    // Swap back
    X_input.swap(X_system);
    R_input.swap(R_system);

    // No errors, we hope
    return 0;
  }


  PetscErrorCode
  __libmesh_petsc_diff_solver_jacobian (SNES,
                                        Vec x,
#if PETSC_RELEASE_LESS_THAN(3,5,0)
                                        Mat * libmesh_dbg_var(j),
                                        Mat * pc,
                                        MatStructure * msflag,
#else
                                        Mat libmesh_dbg_var(j),
                                        Mat pc,
#endif
                                        void * ctx)
  {
    libmesh_assert(x);
    libmesh_assert(j);
    //  libmesh_assert_equal_to (pc, j);  // We don't use separate preconditioners yet
    libmesh_assert(ctx);

    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver*> (ctx));
    ImplicitSystem & sys = solver.system();

    if (solver.verbose)
      libMesh::out << "Assembling the Jacobian" << std::endl;

    PetscVector<Number> & X_system =
      *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X_input(x, sys.comm());

#if PETSC_RELEASE_LESS_THAN(3,5,0)
    PetscMatrix<Number> J_input(*pc, sys.comm());
#else
    PetscMatrix<Number> J_input(pc, sys.comm());
#endif
    PetscMatrix<Number> & J_system =
      *cast_ptr<PetscMatrix<Number> *>(sys.matrix);

    // DiffSystem assembles from the solution and into the jacobian, so
    // swap those with our input vectors before assembling.  They'll
    // probably already be references to the same vectors, but PETSc
    // might do something tricky.
    X_input.swap(X_system);
    J_input.swap(J_system);

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(false, true);
    J_system.close();

    // Swap back
    X_input.swap(X_system);
    J_input.swap(J_system);

#if PETSC_RELEASE_LESS_THAN(3,5,0)
    *msflag = SAME_NONZERO_PATTERN;
#endif
    // No errors, we hope
    return 0;
  }

} // extern "C"


PetscDiffSolver::PetscDiffSolver (sys_type & s)
  : Parent(s)
{
}


void PetscDiffSolver::init ()
{
  LOG_SCOPE("init()", "PetscDiffSolver");

  Parent::init();

  this->setup_petsc_data();
}



PetscDiffSolver::~PetscDiffSolver ()
{
  this->clear();
}



void PetscDiffSolver::clear()
{
  LOG_SCOPE("clear()", "PetscDiffSolver");

  int ierr = LibMeshSNESDestroy(&_snes);
  LIBMESH_CHKERR(ierr);
}



void PetscDiffSolver::reinit()
{
  LOG_SCOPE("reinit()", "PetscDiffSolver");

  // We need to wipe out all the old PETSc data
  // if we are reinit'ing, since we'll need to build
  // it all back up again.
  this->clear();

  Parent::reinit();

  this->setup_petsc_data();
}



DiffSolver::SolveResult convert_solve_result(SNESConvergedReason r)
{
  switch (r)
    {
    case SNES_CONVERGED_FNORM_ABS:
      return DiffSolver::CONVERGED_ABSOLUTE_RESIDUAL;
    case SNES_CONVERGED_FNORM_RELATIVE:
      return DiffSolver::CONVERGED_RELATIVE_RESIDUAL;
#if PETSC_VERSION_LESS_THAN(3,2,1)
    case SNES_CONVERGED_PNORM_RELATIVE:
#else
    case SNES_CONVERGED_SNORM_RELATIVE:
#endif
      return DiffSolver::CONVERGED_RELATIVE_STEP;
    case SNES_CONVERGED_ITS:
    case SNES_CONVERGED_TR_DELTA:
      return DiffSolver::CONVERGED_NO_REASON;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_FNORM_NAN:
#if !PETSC_VERSION_LESS_THAN(3,3,0)
    case SNES_DIVERGED_INNER:
#endif
    case SNES_DIVERGED_LINEAR_SOLVE:
    case SNES_DIVERGED_LOCAL_MIN:
      return DiffSolver::DIVERGED_NO_REASON;
    case SNES_DIVERGED_MAX_IT:
      return DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
#if PETSC_VERSION_LESS_THAN(3,2,0)
    case SNES_DIVERGED_LS_FAILURE:
#else
    case SNES_DIVERGED_LINE_SEARCH:
#endif
      return DiffSolver::DIVERGED_BACKTRACKING_FAILURE;
      // In PETSc, SNES_CONVERGED_ITERATING means
      // the solve is still iterating, but by the
      // time we get here, we must have either
      // converged or diverged, so
      // SNES_CONVERGED_ITERATING is invalid.
    case SNES_CONVERGED_ITERATING:
      return DiffSolver::INVALID_SOLVE_RESULT;
    default:
      break;
    }
  return DiffSolver::INVALID_SOLVE_RESULT;
}



unsigned int PetscDiffSolver::solve()
{
  LOG_SCOPE("solve()", "PetscDiffSolver");

  PetscVector<Number> & x =
    *(cast_ptr<PetscVector<Number> *>(_system.solution.get()));
  PetscMatrix<Number> & jac =
    *(cast_ptr<PetscMatrix<Number> *>(_system.matrix));
  PetscVector<Number> & r =
    *(cast_ptr<PetscVector<Number> *>(_system.rhs));

  int ierr = 0;

  ierr = SNESSetFunction (_snes, r.vec(),
                          __libmesh_petsc_diff_solver_residual, this);
  LIBMESH_CHKERR(ierr);

  ierr = SNESSetJacobian (_snes, jac.mat(), jac.mat(),
                          __libmesh_petsc_diff_solver_jacobian, this);
  LIBMESH_CHKERR(ierr);

  ierr = SNESSolve (_snes, PETSC_NULL, x.vec());
  LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  SNESConvergedReason reason;
  SNESGetConvergedReason(_snes, &reason);

  return convert_solve_result(reason);
}

void PetscDiffSolver::setup_petsc_data()
{
  int ierr=0;

  ierr = SNESCreate(this->comm().get(),&_snes);
  LIBMESH_CHKERR(ierr);

  ierr = SNESMonitorSet (_snes, __libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
  LIBMESH_CHKERR(ierr);

  if (libMesh::on_command_line("--solver-system-names"))
    {
      ierr = SNESSetOptionsPrefix(_snes, (_system.name()+"_").c_str());
      LIBMESH_CHKERR(ierr);
    }

  ierr = SNESSetFromOptions(_snes);
  LIBMESH_CHKERR(ierr);

  KSP my_ksp;
  ierr = SNESGetKSP(_snes, &my_ksp);
  LIBMESH_CHKERR(ierr);

  PC my_pc;
  ierr = KSPGetPC(my_ksp, &my_pc);
  LIBMESH_CHKERR(ierr);

  petsc_auto_fieldsplit(my_pc, _system);
}

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
