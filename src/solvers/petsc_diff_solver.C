// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
  _libmesh_petsc_diff_solver_monitor (SNES snes,
                                      PetscInt its,
                                      PetscReal fnorm,
                                      void * ctx)
  {
    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver *> (ctx));

    if (solver.verbose)
      libMesh::out << "  PetscDiffSolver step " << its
                   << ", |residual|_2 = " << fnorm << std::endl;
    if (solver.linear_solution_monitor.get())
    {
      PetscErrorCode ierr = 0;

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
      ierr = SNESGetFunction(snes, &petsc_res, nullptr, nullptr);
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
  _libmesh_petsc_diff_solver_residual (SNES, Vec x, Vec r, void * ctx)
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
    // something tricky. ... or we might do something tricky.  If
    // we're solving only on constrained DoFs, we'll need to do some
    // scatters here.
    if (solver.get_restrict_solves_to_unconstrained())
      {
        VecScatterBeginEnd(solver.comm(), solver.scatter,
                           X_input.vec(), X_system.vec(),
                           INSERT_VALUES, SCATTER_REVERSE);
      }
    else
      {
        X_input.swap(X_system);
        R_input.swap(R_system);
      }

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(true, false);
    R_system.close();

    // Swap back
    if (solver.get_restrict_solves_to_unconstrained())
      {
        VecScatterBeginEnd(solver.comm(), solver.scatter,
                           R_system.vec(), R_input.vec(),
                           INSERT_VALUES, SCATTER_FORWARD);
      }
    else
      {
        X_input.swap(X_system);
        R_input.swap(R_system);
      }

    // No errors, we hope
    return 0;
  }


  PetscErrorCode
  _libmesh_petsc_diff_solver_jacobian (SNES,
                                       Vec x,
                                       Mat libmesh_dbg_var(j),
                                       Mat pc,
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

    PetscMatrix<Number> J_input(pc, sys.comm());
    PetscMatrix<Number> & J_system =
      *cast_ptr<PetscMatrix<Number> *>(sys.matrix);

    // DiffSystem assembles from the solution and into the jacobian, so
    // swap those with our input vectors before assembling.  They'll
    // probably already be references to the same vectors, but PETSc
    // might do something tricky.  ... or we might do something
    // tricky.  If we're solving only on constrained DoFs, we'll need
    // to do some scatters here.
    if (solver.get_restrict_solves_to_unconstrained())
      {
        VecScatterBeginEnd(solver.comm(), solver.scatter,
                           X_input.vec(), X_system.vec(),
                           INSERT_VALUES, SCATTER_REVERSE);
      }
    else
      {
        X_input.swap(X_system);
        J_input.swap(J_system);
      }

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(false, true);
    J_system.close();

    // Swap back
    if (solver.get_restrict_solves_to_unconstrained())
      {
        PetscInt ierr =
          LibMeshCreateSubMatrix(J_system.mat(),
                                 solver._unconstrained_dofs_is,
                                 solver._unconstrained_dofs_is,
                                 solver.submat_created ?
                                 MAT_REUSE_MATRIX :
                                 MAT_INITIAL_MATRIX,
                                 solver.submat.get());
        LIBMESH_CHKERR(ierr);
        solver.submat_created = true;
      }
    else
      {
        X_input.swap(X_system);
        J_input.swap(J_system);
      }

    // No errors, we hope
    return 0;
  }

} // extern "C"


PetscDiffSolver::PetscDiffSolver (sys_type & s)
  : Parent(s), submat_created(false),
  _restrict_to_unconstrained(false)
{
}


void PetscDiffSolver::init ()
{
  LOG_SCOPE("init()", "PetscDiffSolver");

  Parent::init();

  this->setup_petsc_data();
}



PetscDiffSolver::~PetscDiffSolver () = default;



void PetscDiffSolver::clear()
{
  LOG_SCOPE("clear()", "PetscDiffSolver");

  // calls custom deleter
  _snes.destroy();

#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
  _dm_wrapper.clear();
#endif
#endif
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
    case SNES_CONVERGED_SNORM_RELATIVE:
      return DiffSolver::CONVERGED_RELATIVE_STEP;
    case SNES_CONVERGED_ITS:
      // SNES_CONVERGED_TR_DELTA was changed to a diverged condition,
      // SNES_DIVERGED_TR_DELTA, in PETSc 1c6b2ff8df. This change will
      // likely be in 3.12 and later releases.
#if PETSC_RELEASE_LESS_THAN(3,12,0)
    case SNES_CONVERGED_TR_DELTA:
#endif
      return DiffSolver::CONVERGED_NO_REASON;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_FNORM_NAN:
    case SNES_DIVERGED_INNER:
    case SNES_DIVERGED_LINEAR_SOLVE:
    case SNES_DIVERGED_LOCAL_MIN:
      return DiffSolver::DIVERGED_NO_REASON;
    case SNES_DIVERGED_MAX_IT:
      return DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
    case SNES_DIVERGED_LINE_SEARCH:
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


void PetscDiffSolver::restrict_solves_to_unconstrained(bool restricting)
{
  _restrict_to_unconstrained = restricting;
}


bool PetscDiffSolver::get_restrict_solves_to_unconstrained()
{
  return _restrict_to_unconstrained;
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

  WrappedPetsc<Vec> subrhs, subsolution;

  PetscErrorCode ierr = 0;

  if (_restrict_to_unconstrained)
    {
      std::vector<PetscInt> unconstrained_dofs;
      const DofMap & dofmap = this->system().get_dof_map();
      for (dof_id_type i = dofmap.first_dof(),
                   end_i = dofmap.end_dof();
           i != end_i; ++i)
        if (!dofmap.is_constrained_dof(i))
          unconstrained_dofs.push_back(cast_int<PetscInt>(i));

      const PetscInt is_local_size =
        cast_int<PetscInt>(unconstrained_dofs.size());

      ierr = ISCreateGeneral(this->comm().get(),
                             cast_int<PetscInt>(unconstrained_dofs.size()),
                             unconstrained_dofs.data(), PETSC_COPY_VALUES,
                             _unconstrained_dofs_is.get());

      ierr = VecCreate(this->comm().get(), subrhs.get());
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subrhs, is_local_size, PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subrhs);
      LIBMESH_CHKERR(ierr);

      ierr = VecCreate(this->comm().get(), subsolution.get());
      LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(subsolution, is_local_size, PETSC_DECIDE);
      LIBMESH_CHKERR(ierr);
      ierr = VecSetFromOptions(subsolution);
      LIBMESH_CHKERR(ierr);

      ierr = VecScatterCreate(r.vec(), _unconstrained_dofs_is,
                              subrhs, nullptr, scatter.get());
      LIBMESH_CHKERR(ierr);

      VecScatterBeginEnd(this->comm(), scatter, r.vec(), subrhs, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterBeginEnd(this->comm(), scatter, x.vec(), subsolution, INSERT_VALUES, SCATTER_FORWARD);

      LIBMESH_CHKERR(ierr);

      ierr = SNESSetFunction (_snes, subrhs,
                              _libmesh_petsc_diff_solver_residual, this);
      LIBMESH_CHKERR(ierr);

      ierr = SNESSetJacobian (_snes, submat, submat,
                              _libmesh_petsc_diff_solver_jacobian, this);
    }
  else
    {
      ierr = SNESSetFunction (_snes, r.vec(),
                              _libmesh_petsc_diff_solver_residual, this);
      LIBMESH_CHKERR(ierr);

      ierr = SNESSetJacobian (_snes, jac.mat(), jac.mat(),
                              _libmesh_petsc_diff_solver_jacobian, this);
    }
  LIBMESH_CHKERR(ierr);

  ierr = SNESSetFromOptions(_snes);
  LIBMESH_CHKERR(ierr);

  if (_restrict_to_unconstrained)
    ierr = SNESSolve (_snes, PETSC_NULL, x.vec());
  else
    ierr = SNESSolve (_snes, PETSC_NULL, x.vec());
  LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  SNESConvergedReason reason;
  SNESGetConvergedReason(_snes, &reason);

  PetscInt l_its, nl_its;
  ierr = SNESGetLinearSolveIterations(_snes, &l_its);
  LIBMESH_CHKERR(ierr);
  this->_inner_iterations = l_its;

  ierr = SNESGetIterationNumber(_snes, &nl_its);
  LIBMESH_CHKERR(ierr);
  this->_outer_iterations = nl_its;

  return convert_solve_result(reason);
}

void PetscDiffSolver::setup_petsc_data()
{
  PetscErrorCode ierr = 0;

  ierr = SNESCreate(this->comm().get(), _snes.get());
  LIBMESH_CHKERR(ierr);

  ierr = SNESMonitorSet (_snes, _libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
  LIBMESH_CHKERR(ierr);

  if (libMesh::on_command_line("--solver-system-names"))
    {
      ierr = SNESSetOptionsPrefix(_snes, (_system.name()+"_").c_str());
      LIBMESH_CHKERR(ierr);
    }

  bool use_petsc_dm = libMesh::on_command_line("--use_petsc_dm");

  // This needs to be called before SNESSetFromOptions
#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
  if (use_petsc_dm)
    this->_dm_wrapper.init_and_attach_petscdm(_system, *_snes);
#endif
#endif

  // If we're not using PETSc DM, let's keep around
  // the old style for fieldsplit
  if (!use_petsc_dm)
    {
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
}

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
