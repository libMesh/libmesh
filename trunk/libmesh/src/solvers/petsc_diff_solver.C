
#include "diff_system.h"
#include "dof_map.h"
#include "libmesh_logging.h"
#include "petsc_diff_solver.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"

#ifdef HAVE_PETSC

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
__libmesh_petsc_diff_solver_monitor (SNES, PetscInt its,
                                     PetscReal fnorm, void *)
{
  std::cout << "  PetscDiffSolver step " << its
            << ", |residual|_2 = " << fnorm << std::endl;

  return 0;
}

// Functions to hand to PETSc's SNES,
// which compute the residual or jacobian at X
PetscErrorCode
__libmesh_petsc_diff_solver_residual (SNES, Vec x, Vec r, void *ctx)
{
  assert (x   != NULL);
  assert (r   != NULL);
  assert (ctx != NULL);

  PetscDiffSolver& solver =
    *(static_cast<PetscDiffSolver*> (ctx));
  DifferentiableSystem &sys = solver.system();

  PetscVector<Number> X_input(x), R_input(r);
  PetscVector<Number>& X_system =
    *dynamic_cast<PetscVector<Number>*>(sys.solution.get());
  PetscVector<Number>& R_system =
    *dynamic_cast<PetscVector<Number>*>(sys.rhs);

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

  // Swap back
  X_input.swap(X_system);
  R_input.swap(R_system);

  // No errors, we hope
  return 0;
}


PetscErrorCode
__libmesh_petsc_diff_solver_jacobian (SNES, Vec x, Mat *j, Mat *pc,
                                      MatStructure *msflag, void *ctx)
{
  assert (x   != NULL);
  assert (j   != NULL);
//  assert (pc  == j);  // We don't use separate preconditioners yet
  assert (ctx != NULL);

  PetscDiffSolver& solver =
    *(static_cast<PetscDiffSolver*> (ctx));
  DifferentiableSystem &sys = solver.system();

  PetscVector<Number> X_input(x);
  PetscVector<Number>& X_system =
    *dynamic_cast<PetscVector<Number>*>(sys.solution.get());

  PetscMatrix<Number> J_input(*pc);
  PetscMatrix<Number>& J_system =
    *dynamic_cast<PetscMatrix<Number>*>(sys.matrix);

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
  int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,1,2)
  // At least until Petsc 2.1.1, the SNESCreate had a different
  // calling syntax.  The second argument was of type SNESProblemType,
  // and could have a value of either SNES_NONLINEAR_EQUATIONS or
  // SNES_UNCONSTRAINED_MINIMIZATION.
  ierr = SNESCreate(libMesh::COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &_snes);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  ierr = SNESCreate(libMesh::COMM_WORLD,&_snes);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = SNESSetMonitor (_snes, __libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
#else
  // API name change in PETSc 2.3.3
  ierr = SNESMonitorSet (_snes, __libmesh_petsc_diff_solver_monitor,
                         this, PETSC_NULL);
#endif
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SNESSetFromOptions(_snes);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

}



PetscDiffSolver::~PetscDiffSolver ()
{
  int ierr=0;

  ierr = SNESDestroy(_snes);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
}



void PetscDiffSolver::reinit()
{
  Parent::reinit();
}



unsigned int PetscDiffSolver::solve()
{
  START_LOG("solve()", "PetscDiffSolver");

  PetscVector<Number> &x =
    *(dynamic_cast<PetscVector<Number>*>(_system.solution.get()));
  PetscMatrix<Number> &jac =
    *(dynamic_cast<PetscMatrix<Number>*>(_system.matrix));
  PetscVector<Number> &r =
    *(dynamic_cast<PetscVector<Number>*>(_system.rhs));

  x.close();
  r.close();
  jac.close();

#ifdef ENABLE_AMR
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  int ierr = 0;

  ierr = SNESSetFunction (_snes, r.vec(),
                          __libmesh_petsc_diff_solver_residual, this);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SNESSetJacobian (_snes, jac.mat(), jac.mat(),
                          __libmesh_petsc_diff_solver_jacobian, this);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

# if PETSC_VERSION_LESS_THAN(2,2,0)

  ierr = SNESSolve (_snes, x.vec(), &_outer_iterations);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.x style
#elif PETSC_VERSION_LESS_THAN(2,3,0)

  ierr = SNESSolve (_snes, x.vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.3.x & newer style
#else

  ierr = SNESSolve (_snes, PETSC_NULL, x.vec());
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif

  STOP_LOG("solve()", "PetscDiffSolver");

  // FIXME - We'll worry about getting the solve result right later...
  
  return DiffSolver::CONVERGED_RELATIVE_RESIDUAL;
}

#endif // HAVE_PETSC
