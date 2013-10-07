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


#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_diff_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

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
  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

// Function to hand to PETSc's SNES,
// which monitors convergence at X
PetscErrorCode
__libmesh_petsc_diff_solver_monitor (SNES snes, PetscInt its,
                                     PetscReal fnorm, void *ctx)
{
  PetscDiffSolver& solver =
    *(static_cast<PetscDiffSolver*> (ctx));

  if (solver.verbose)
    libMesh::out << "  PetscDiffSolver step " << its
                 << ", |residual|_2 = " << fnorm << std::endl;
  if (solver.linear_solution_monitor.get()) {
    int ierr = 0;

    Vec petsc_delta_u;
    ierr = SNESGetSolutionUpdate(snes, &petsc_delta_u);
    CHKERRABORT(libMesh::COMM_WORLD, ierr);
    PetscVector<Number> delta_u(petsc_delta_u);
    delta_u.close();

    Vec petsc_u;
    ierr = SNESGetSolution(snes, &petsc_u);
    CHKERRABORT(libMesh::COMM_WORLD, ierr);
    PetscVector<Number> u(petsc_u);
    u.close();

    Vec petsc_res;
    ierr = SNESGetFunction(snes, &petsc_res, NULL, NULL);
    CHKERRABORT(libMesh::COMM_WORLD, ierr);
    PetscVector<Number> res(petsc_res);
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
__libmesh_petsc_diff_solver_residual (SNES, Vec x, Vec r, void *ctx)
{
  libmesh_assert(x);
  libmesh_assert(r);
  libmesh_assert(ctx);

  PetscDiffSolver& solver =
    *(static_cast<PetscDiffSolver*> (ctx));
  ImplicitSystem &sys = solver.system();

  if (solver.verbose)
    libMesh::out << "Assembling the residual" << std::endl;

  PetscVector<Number>& X_system =
    *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
  PetscVector<Number>& R_system =
    *libmesh_cast_ptr<PetscVector<Number>*>(sys.rhs);
  PetscVector<Number> X_input(x, sys.comm()), R_input(r, sys.comm());

  // DiffSystem assembles from the solution and into the rhs, so swap
  // those with our input vectors before assembling.  They'll probably
  // already be references to the same vectors, but PETSc might do
  // something tricky.
  X_input.swap(X_system);
  R_input.swap(R_system);

  // We may need to correct a non-conforming solution
  sys.get_dof_map().enforce_constraints_exactly(sys);

  // We may need to localize a parallel solution
  sys.update();

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
__libmesh_petsc_diff_solver_jacobian (SNES, Vec x, Mat *libmesh_dbg_var(j), Mat *pc,
                                      MatStructure *msflag, void *ctx)
{
  libmesh_assert(x);
  libmesh_assert(j);
//  libmesh_assert_equal_to (pc, j);  // We don't use separate preconditioners yet
  libmesh_assert(ctx);

  PetscDiffSolver& solver =
    *(static_cast<PetscDiffSolver*> (ctx));
  ImplicitSystem &sys = solver.system();

  if (solver.verbose)
    libMesh::out << "Assembling the Jacobian" << std::endl;

  PetscVector<Number>& X_system =
    *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
  PetscVector<Number> X_input(x, sys.comm());

  PetscMatrix<Number> J_input(*pc, sys.comm());
  PetscMatrix<Number>& J_system =
    *libmesh_cast_ptr<PetscMatrix<Number>*>(sys.matrix);

  // DiffSystem assembles from the solution and into the jacobian, so
  // swap those with our input vectors before assembling.  They'll
  // probably already be references to the same vectors, but PETSc
  // might do something tricky.
  X_input.swap(X_system);
  J_input.swap(J_system);

  // We may need to correct a non-conforming solution
  sys.get_dof_map().enforce_constraints_exactly(sys);

  // We may need to localize a parallel solution
  sys.update();

  // Do DiffSystem assembly
  sys.assembly(false, true);
  J_system.close();

  // Swap back
  X_input.swap(X_system);
  J_input.swap(J_system);

  *msflag = SAME_NONZERO_PATTERN;

  // No errors, we hope
  return 0;
}

} // extern "C"


PetscDiffSolver::PetscDiffSolver (sys_type& s)
  : Parent(s)
{
}


void PetscDiffSolver::init ()
{
  START_LOG("init()", "PetscDiffSolver");

  Parent::init();

  int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,1,2)
  // At least until Petsc 2.1.1, the SNESCreate had a different
  // calling syntax.  The second argument was of type SNESProblemType,
  // and could have a value of either SNES_NONLINEAR_EQUATIONS or
  // SNES_UNCONSTRAINED_MINIMIZATION.
  ierr = SNESCreate(this->comm().get(), SNES_NONLINEAR_EQUATIONS, &_snes);
  LIBMESH_CHKERRABORT(ierr);
#else
  ierr = SNESCreate(this->comm().get(),&_snes);
  LIBMESH_CHKERRABORT(ierr);
#endif

#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = SNESSetMonitor (_snes, __libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
#else
  // API name change in PETSc 2.3.3
  ierr = SNESMonitorSet (_snes, __libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
#endif
  LIBMESH_CHKERRABORT(ierr);

  ierr = SNESSetFromOptions(_snes);
  LIBMESH_CHKERRABORT(ierr);

  STOP_LOG("init()", "PetscDiffSolver");
}



PetscDiffSolver::~PetscDiffSolver ()
{
}



void PetscDiffSolver::clear()
{
  START_LOG("clear()", "PetscDiffSolver");

  int ierr=0;

  ierr = LibMeshSNESDestroy(&_snes);
  LIBMESH_CHKERRABORT(ierr);

  STOP_LOG("clear()", "PetscDiffSolver");
}



void PetscDiffSolver::reinit()
{
  Parent::reinit();
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
#if !PETSC_VERSION_LESS_THAN(2,3,3)
    case SNES_CONVERGED_ITS:
#endif
    case SNES_CONVERGED_TR_DELTA:
      return DiffSolver::CONVERGED_NO_REASON;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_FNORM_NAN:
#if !PETSC_VERSION_LESS_THAN(3,3,0)
    case SNES_DIVERGED_INNER:
#endif
#if !PETSC_VERSION_LESS_THAN(2,3,2)
    case SNES_DIVERGED_LINEAR_SOLVE:
#endif
    case SNES_DIVERGED_LOCAL_MIN:
      return DiffSolver::DIVERGED_NO_REASON;
    case SNES_DIVERGED_MAX_IT:
      return DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
#if !PETSC_VERSION_LESS_THAN(3,2,0)
    case SNES_DIVERGED_LINE_SEARCH:
      return DiffSolver::DIVERGED_BACKTRACKING_FAILURE;
#endif
      // In PETSc, SNES_CONVERGED_ITERATING means
      // the solve is still iterating, but by the
      // time we get here, we must have either
      // converged or diverged, so
      // SNES_CONVERGED_ITERATING is invalid.
    case SNES_CONVERGED_ITERATING:
      return DiffSolver::INVALID_SOLVE_RESULT;
    }
  return DiffSolver::INVALID_SOLVE_RESULT;
}



unsigned int PetscDiffSolver::solve()
{
  this->init();

  START_LOG("solve()", "PetscDiffSolver");

  PetscVector<Number> &x =
    *(libmesh_cast_ptr<PetscVector<Number>*>(_system.solution.get()));
  PetscMatrix<Number> &jac =
    *(libmesh_cast_ptr<PetscMatrix<Number>*>(_system.matrix));
  PetscVector<Number> &r =
    *(libmesh_cast_ptr<PetscVector<Number>*>(_system.rhs));

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  int ierr = 0;

  ierr = SNESSetFunction (_snes, r.vec(),
                          __libmesh_petsc_diff_solver_residual, this);
    LIBMESH_CHKERRABORT(ierr);

  ierr = SNESSetJacobian (_snes, jac.mat(), jac.mat(),
                          __libmesh_petsc_diff_solver_jacobian, this);
    LIBMESH_CHKERRABORT(ierr);

# if PETSC_VERSION_LESS_THAN(2,2,0)

  ierr = SNESSolve (_snes, x.vec(), &_outer_iterations);
         LIBMESH_CHKERRABORT(ierr);

// 2.2.x style
#elif PETSC_VERSION_LESS_THAN(2,3,0)

  ierr = SNESSolve (_snes, x.vec());
         LIBMESH_CHKERRABORT(ierr);

// 2.3.x & newer style
#else

  ierr = SNESSolve (_snes, PETSC_NULL, x.vec());
         LIBMESH_CHKERRABORT(ierr);

#endif

  STOP_LOG("solve()", "PetscDiffSolver");

  SNESConvergedReason reason;
  SNESGetConvergedReason(_snes, &reason);

  this->clear();

  return convert_solve_result(reason);
}


} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
