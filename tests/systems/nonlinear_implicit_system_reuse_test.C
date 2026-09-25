// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include <libmesh/libmesh_config.h>
#include <libmesh/petsc_macro.h> // defines PETSC_VERSION_LESS_THAN() without PETSc too

// The code under test reads the SNES with SNESGetUseMatrixFree(), which exists from PETSc
// 3.8 on; below that, and without PETSc, the fix compiles out and so does this file.
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,8,0)

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/nonlinear_solver.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/parallel_implementation.h> // max()
#include <libmesh/petsc_nonlinear_solver.h>
#include <libmesh/solver_configuration.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/wrapped_petsc.h>
#include <libmesh/int_range.h>

#include <cmath>
#include <string>
#include <vector>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
// R_i(x) = x_i^3 + x_i - 1, a diagonal, everywhere nonsingular nonlinear system:
// J_ii = 3 x_i^2 + 1 >= 1, so it converges from x = 0 in a handful of Newton steps
// and needs no assembly machinery. Its single real root, to 17 digits:
const Real cubic_root = 0.68232780382801932;

// The residual is an object, not a function pointer: the matrix-free residual
// callback (libmesh_petsc_snes_mffd_residual) serves mffd_residual_object, then
// residual_object, and never the residual function pointer.
class CubicResidual : public NonlinearImplicitSystem::ComputeResidual
{
public:
  virtual void residual(const NumericVector<Number> & X,
                        NumericVector<Number> & R,
                        NonlinearImplicitSystem &) override
  {
    for (auto i : make_range(X.first_local_index(), X.last_local_index()))
      R.set(i, X(i) * X(i) * X(i) + X(i) - Number(1));
    R.close();
  }
};

class CubicJacobian : public NonlinearImplicitSystem::ComputeJacobian
{
public:
  virtual void jacobian(const NumericVector<Number> & X,
                        SparseMatrix<Number> & J,
                        NonlinearImplicitSystem &) override
  {
    J.zero();
    for (auto i : make_range(X.first_local_index(), X.last_local_index()))
      J.set(i, i, Number(3) * X(i) * X(i) + Number(1));
    J.close();
  }
};

// R_i(x) = x_i, for the left nonlinear preconditioner case: zero at x = 0, so a solve
// started there converges at its initial iterate, before any Jacobian is needed.
class IdentityResidual : public NonlinearImplicitSystem::ComputeResidual
{
public:
  virtual void residual(const NumericVector<Number> & X,
                        NumericVector<Number> & R,
                        NonlinearImplicitSystem &) override
  {
    for (auto i : make_range(X.first_local_index(), X.last_local_index()))
      R.set(i, X(i));
    R.close();
  }
};

// Sets or clears PETSc options for one case and, when it goes out of scope, puts the
// options database back the way it was. The driver runs every suite in one process, so an
// option left behind would reach the cases and suites that follow (a case that fails
// throws past any cleanup written after it), and an option that was already there
// when the case began must come back, not disappear. Every option is recorded before
// the first change, so the guard is declared before the case sets anything.
class PetscOptionsScope
{
public:
  explicit PetscOptionsScope(const std::vector<std::string> & names)
  {
    for (const auto & name : names)
      {
        SavedOption saved;
        saved.name = name;
        char value[PETSC_MAX_PATH_LEN];
        PetscBool was_set = PETSC_FALSE;
        LibmeshPetscCallA(TestCommWorld->get(),
                          PetscOptionsGetString(LIBMESH_PETSC_NULLPTR, LIBMESH_PETSC_NULLPTR,
                                                name.c_str(), value, sizeof(value), &was_set));
        saved.was_set = (was_set == PETSC_TRUE);
        if (saved.was_set)
          saved.value = value; // empty for an option that was given without a value
        _saved.push_back(saved);
      }
  }

  ~PetscOptionsScope()
  {
    // A destructor must not throw, so the PETSc return codes are ignored here. An
    // empty value puts back an option that was given without one.
    for (const auto & saved : _saved)
      {
        if (saved.was_set)
          libmesh_ignore(PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, saved.name.c_str(),
                                              saved.value.c_str()));
        else
          libmesh_ignore(PetscOptionsClearValue(LIBMESH_PETSC_NULLPTR, saved.name.c_str()));
      }
  }

  // Sets an option that was named to the constructor.
  void set(const std::string & name, const std::string & value)
  {
    check_recorded(name);
    LibmeshPetscCallA(TestCommWorld->get(),
                      PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, name.c_str(), value.c_str()));
  }

  // Removes an option that was named to the constructor, so that a value left in the
  // database by whatever ran before cannot reach the case; scope exit puts it back.
  void clear(const std::string & name)
  {
    check_recorded(name);
    LibmeshPetscCallA(TestCommWorld->get(),
                      PetscOptionsClearValue(LIBMESH_PETSC_NULLPTR, name.c_str()));
  }

private:
  void check_recorded(const std::string & name) const
  {
    bool recorded = false;
    for (const auto & saved : _saved)
      if (saved.name == name)
        recorded = true;
    libmesh_error_msg_if(!recorded,
                         "PetscOptionsScope: " << name << " was not named to the constructor");
  }

  struct SavedOption
  {
    std::string name;
    bool was_set = false;
    std::string value;
  };

  std::vector<SavedOption> _saved;
};

// For the replaced-preconditioner case. The hook installs a brand new shell PC on
// every solve and, on the second solve, a brand new KSP as well: KSPSetPC() releases
// the KSP's reference to the previous PC and SNESSetKSP() the SNES's reference to
// the previous KSP, so a KSP or PC pointer the solver read before the hook ran is
// stale after it. The first solve's KSP is kept alive here, together with the PC it
// holds, so that a solver acting through a stale pointer changes an object the test
// can inspect instead of freed memory, about which an optimized build says nothing.
// An identity apply is enough: the system is diagonal.
PetscErrorCode shell_pc_apply(PC, Vec x, Vec y) { return VecCopy(x, y); }

class ReplaceKSPAndPCConfiguration : public SolverConfiguration
{
public:
  explicit ReplaceKSPAndPCConfiguration(PetscNonlinearSolver<Number> & solver) : _solver(solver) {}

  virtual void configure_solver() override
  {
    const Parallel::Communicator & comm = _solver.comm();
    // solve() has already called init(), so snes() is a plain accessor here.
    SNES snes = _solver.snes();

    if (n_configures == 1)
      {
        // Second solve: keep the first solve's KSP alive, then replace it. The new KSP
        // gets its operators from the Newton iteration and its PC from the code below.
        LibmeshPetscCall2(comm, SNESGetKSP(snes, previous_ksp.get()));
        LibmeshPetscCall2(comm, PetscObjectReference((PetscObject)(*previous_ksp)));

        WrappedPetsc<KSP> new_ksp;
        LibmeshPetscCall2(comm, KSPCreate(comm.get(), new_ksp.get()));
        LibmeshPetscCall2(comm, SNESSetKSP(snes, new_ksp)); // the SNES takes its own reference
      }

    KSP ksp;
    WrappedPetsc<PC> pc;
    LibmeshPetscCall2(comm, SNESGetKSP(snes, &ksp));
    LibmeshPetscCall2(comm, PCCreate(comm.get(), pc.get()));
    LibmeshPetscCall2(comm, PCSetType(pc, PCSHELL));
    LibmeshPetscCall2(comm, PCShellSetApply(pc, shell_pc_apply));
    LibmeshPetscCall2(comm, KSPSetPC(ksp, pc)); // the KSP takes its own reference
    // installed_pc is borrowed, for the assertions only.
    LibmeshPetscCall2(comm, KSPGetPC(ksp, &installed_pc));
    ++n_configures;
  }

  WrappedPetsc<KSP> previous_ksp;
  PC installed_pc = LIBMESH_PETSC_NULLPTR;
  unsigned int n_configures = 0;

private:
  PetscNonlinearSolver<Number> & _solver;
};
}

class NonlinearImplicitSystemReuseTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( NonlinearImplicitSystemReuseTest );
  CPPUNIT_TEST( testReuseKeepsMatrixFreeOperator );
  CPPUNIT_TEST( testReuseKeepsFullyMatrixFreeOperator );
  CPPUNIT_TEST( testReuseResidualOnlyMatrixFree );
  CPPUNIT_TEST( testReuseWithReplacedPreconditioner );
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  CPPUNIT_TEST( testReuseRejectsModeChange );
#endif
  CPPUNIT_TEST( testReuseKeepsMatrixFreeUnderFiniteDifferenceOption );
#if defined(LIBMESH_ENABLE_EXCEPTIONS) && !defined(PETSC_USE_DEBUG)
  CPPUNIT_TEST( testReuseWithLeftNonlinearPreconditioner );
#endif
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}

  // Nine EDGE2 elements with one unknown each: the problem is diagonal and needs nine
  // unknowns, not geometry, so the cases run on a libMesh built for one dimension.
  // Keep the SNES (and so its matrix-free operator) alive across solves: this has to go
  // through the EquationSystems parameters, because NonlinearImplicitSystem::solve()
  // copies them onto the solver on every solve. The iteration limit is large enough
  // that libMesh's recalculation monitor never arms a rebuild; the rebuild policy is
  // not what these cases are about.
  static NonlinearImplicitSystem & build(EquationSystems & es, UnstructuredMesh & mesh)
  {
    MeshTools::Generation::build_line(mesh, 9, 0., 1., EDGE2);
    NonlinearImplicitSystem & system = es.add_system<NonlinearImplicitSystem>("nl");
    system.add_variable("u", CONSTANT, MONOMIAL);
    es.parameters.set<bool>("reuse preconditioner") = true;
    es.parameters.set<unsigned int>
      ("reuse preconditioner maximum linear iterations") = 1000;
    es.init();
    return system;
  }

  // The verdict of an assertion, made the same on every rank: an assertion that failed
  // on one rank alone would unwind that rank alone, and the ranks would part at the
  // next collective call. Every CPPUNIT_ASSERT below goes through this.
  static bool passes_everywhere(bool passes_here)
  {
    unsigned int failed = passes_here ? 0 : 1;
    TestCommWorld->max(failed);
    return failed == 0;
  }

  // The two matrix-free modes, each with both options explicit, so that an option left
  // in the database by whatever ran before cannot turn one mode into the other. PETSc
  // keeps two flags, mf_operator and mf: -snes_mf_operator implies mf, and an explicit
  // -snes_mf false beside -snes_mf_operator true leaves mf false, which SNESSetUp()
  // takes as no matrix-free setup at all. So a -snes_mf_operator solve sets true/true
  // and a -snes_mf solve false/true.
  static void set_matrix_free_options(PetscOptionsScope & options, bool mf_operator)
  {
    options.set("-snes_mf_operator", mf_operator ? "true" : "false");
    options.set("-snes_mf", "true");
  }

  static PetscNonlinearSolver<Number> & petsc_solver(NonlinearImplicitSystem & system)
  {
    return cast_ref<PetscNonlinearSolver<Number> &>(*system.nonlinear_solver);
  }

  // NOTE: snes() calls init(), which is not a no-op after a reuse solve. It attaches a
  // fresh DM and leaves both operators alone, which is all these assertions need.
  static void get_matrices(PetscNonlinearSolver<Number> & solver, Mat & Amat, Mat & Pmat)
  {
    Amat = LIBMESH_PETSC_NULLPTR;
    Pmat = LIBMESH_PETSC_NULLPTR;
    LibmeshPetscCallA(TestCommWorld->get(),
                      SNESGetJacobian(solver.snes(), &Amat, &Pmat,
                                      LIBMESH_PETSC_NULLPTR, LIBMESH_PETSC_NULLPTR));
  }

  static PC current_pc(PetscNonlinearSolver<Number> & solver)
  {
    KSP ksp;
    PC pc;
    LibmeshPetscCallA(TestCommWorld->get(), SNESGetKSP(solver.snes(), &ksp));
    LibmeshPetscCallA(TestCommWorld->get(), KSPGetPC(ksp, &pc));
    return pc;
  }

  // Whether the solution is the root of x^3 + x = 1, on every rank. 1e-6, not 1e-8: a
  // differenced Jacobian at PETSc's default tolerances lands near 1e-8 of the root, so a
  // tighter bound would test the tolerances, not the fix. The difference is taken in
  // Number, so that an imaginary part counts under complex scalars. Each rank checks the
  // entries it owns; the verdict is reduced.
  static bool solution_is_cubic_root(NonlinearImplicitSystem & system)
  {
    bool solution_is_right_here = true;
    for (auto i : make_range(system.solution->first_local_index(),
                             system.solution->last_local_index()))
      {
        const Real error = std::abs((*system.solution)(i) - Number(cubic_root));
        if (libmesh_isnan(error) || libmesh_isinf(error) || error > 1e-6)
          solution_is_right_here = false;
      }
    return passes_everywhere(solution_is_right_here);
  }

  // Cases 1 and 2.
  // mf_operator == true  -> -snes_mf_operator: Amat matrix-free, Pmat assembled
  // mf_operator == false -> -snes_mf:          Amat and Pmat both matrix-free
  void runTwoSolves(bool mf_operator)
  {
    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf"});
    set_matrix_free_options(options, mf_operator);

    {
      CubicResidual cubic_residual;
      CubicJacobian cubic_jacobian;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &cubic_residual;
      system.nonlinear_solver->jacobian_object = &cubic_jacobian;

      for (unsigned int solve_number = 1; solve_number <= 2; ++solve_number)
        {
          // Start every solve away from the root, so that solve 2 does real Newton work
          // (a solve started at solve 1's converged iterate can return before it ever
          // evaluates the Jacobian, which would leave the callback untested).
          system.solution->zero();
          system.solution->close();
          system.update();
          system.solve();
          CPPUNIT_ASSERT_MESSAGE("the solve did no Newton iteration",
                                 passes_everywhere(system.n_nonlinear_iterations() > 0));

          Mat Amat, Pmat;
          get_matrices(petsc_solver(system), Amat, Pmat);
          CPPUNIT_ASSERT_MESSAGE("the SNES has no operator",
                                 passes_everywhere(Amat != nullptr));
          CPPUNIT_ASSERT_MESSAGE("the SNES has no preconditioning matrix",
                                 passes_everywhere(Pmat != nullptr));

          PetscBool is_mffd = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)Amat, MATMFFD, &is_mffd));
          // Solve 2 is the one that regresses: without the fix the operator here is the
          // assembled preconditioning matrix and a matrix-free solve has silently become
          // a Newton solve on it.
          CPPUNIT_ASSERT_MESSAGE("SNES operator is not matrix-free",
                                 passes_everywhere(is_mffd == PETSC_TRUE));

          if (mf_operator)
            CPPUNIT_ASSERT_MESSAGE("operator and preconditioning matrix are the same object",
                                   passes_everywhere(Amat != Pmat));
          else
            CPPUNIT_ASSERT_MESSAGE("-snes_mf must use one matrix for both",
                                   passes_everywhere(Amat == Pmat));
        }
    }
  }

  // Case 3: residual only, as an application that never assembles a matrix under
  // -snes_mf, and, before PETSc 3.26, an explicit -pc_type that PETSc's matrix-free
  // setup overrides on solve 1 and that nothing overrides on solve 2 without the
  // repair. From 3.26 on the setup leaves a chosen PC type alone (and errors on one
  // that needs an assembled matrix), so the option is not set there and the PCNONE
  // asserted below is the one SNESSetFromOptions() chose on solve 1 and kept. This is
  // the case that asserts a correct answer, not only a matrix type.
  void testReuseResidualOnlyMatrixFree()
  {
    LOG_UNIT_TEST;

    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf", "-pc_type"});
    set_matrix_free_options(options, false);
    // The case starts without a -pc_type: from PETSc 3.26 on an inherited one would be
    // honored and fail the PCNONE assertion below on a correct library.
    options.clear("-pc_type");
#if PETSC_VERSION_LESS_THAN(3,26,0)
    options.set("-pc_type", "ilu");
#endif

    {
      CubicResidual cubic_residual;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &cubic_residual;

      for (unsigned int solve_number = 1; solve_number <= 2; ++solve_number)
        {
          system.solution->zero(); // away from the root: see case 1
          system.solution->close();
          system.update();
          system.solve();
          CPPUNIT_ASSERT_MESSAGE("the solve did no Newton iteration",
                                 passes_everywhere(system.n_nonlinear_iterations() > 0));

          SNESConvergedReason reason;
          LibmeshPetscCallA(TestCommWorld->get(),
                            SNESGetConvergedReason(petsc_solver(system).snes(), &reason));
          CPPUNIT_ASSERT_MESSAGE("the matrix-free solve did not converge",
                                 passes_everywhere(reason > 0));

          CPPUNIT_ASSERT_MESSAGE("the solution is not the root of x^3 + x = 1",
                                 solution_is_cubic_root(system));

          Mat Amat, Pmat;
          get_matrices(petsc_solver(system), Amat, Pmat);
          PetscBool is_mffd = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)Amat, MATMFFD, &is_mffd));
          CPPUNIT_ASSERT_MESSAGE("SNES operator is not matrix-free",
                                 passes_everywhere(is_mffd == PETSC_TRUE));
          CPPUNIT_ASSERT_MESSAGE("-snes_mf must use one matrix for both",
                                 passes_everywhere(Amat == Pmat));

          // -snes_mf runs with no preconditioner. Before PETSc 3.26, without the
          // re-forcing, solve 2 tries to build an ILU from a matrix-free matrix and aborts.
          PetscBool is_none = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)current_pc(petsc_solver(system)),
                                                   PCNONE, &is_none));
          CPPUNIT_ASSERT_MESSAGE("-snes_mf must run with PCNONE",
                                 passes_everywhere(is_none == PETSC_TRUE));
        }
    }
  }

  // Case 4: the PC is replaced on every solve by a configuration hook, and on the
  // second solve the KSP as well, after the solver has already read both once. The
  // repair has to read the KSP and the PC again after the hooks; a pointer read before
  // them is stale. On the second solve those stale pointers are the first solve's KSP
  // and its PC, which the options step has just retyped to ILU (-pc_type is re-read on
  // every solve) and which the hook then unseats. The fixture keeps that KSP alive: a
  // solver acting through a stale pointer would force PCNONE on the PC it still holds,
  // and on freed memory when nothing keeps it alive, which no assertion can see in an
  // optimized build.
  void testReuseWithReplacedPreconditioner()
  {
    LOG_UNIT_TEST;

    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf", "-pc_type"});
    set_matrix_free_options(options, false);
    options.set("-pc_type", "ilu");

    {
      CubicResidual cubic_residual;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &cubic_residual;

      ReplaceKSPAndPCConfiguration configuration(petsc_solver(system));
      system.nonlinear_solver->set_solver_configuration(configuration);

      for (unsigned int solve_number = 1; solve_number <= 2; ++solve_number)
        {
          system.solution->zero(); // away from the root: see case 1
          system.solution->close();
          system.update();
          system.solve();
          CPPUNIT_ASSERT_MESSAGE("the solve did no Newton iteration",
                                 passes_everywhere(system.n_nonlinear_iterations() > 0));

          SNESConvergedReason reason;
          LibmeshPetscCallA(TestCommWorld->get(),
                            SNESGetConvergedReason(petsc_solver(system).snes(), &reason));
          CPPUNIT_ASSERT_MESSAGE("the solve with a replaced preconditioner did not converge",
                                 passes_everywhere(reason > 0));

          // A shell PC is exempt from the PCNONE re-forcing, as it is in PETSc, so the
          // PC of the KSP the SNES holds now must be the one this solve's hook installed.
          PC pc = current_pc(petsc_solver(system));
          PetscBool is_shell = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)pc, PCSHELL, &is_shell));
          CPPUNIT_ASSERT_MESSAGE("the configuration's shell PC did not survive the solve",
                                 passes_everywhere(is_shell == PETSC_TRUE));
          CPPUNIT_ASSERT_MESSAGE("the current KSP's PC is not the one the hook installed",
                                 passes_everywhere(pc == configuration.installed_pc));

          if (solve_number == 2)
            {
              KSP ksp;
              LibmeshPetscCallA(TestCommWorld->get(),
                                SNESGetKSP(petsc_solver(system).snes(), &ksp));
              KSP previous_ksp = configuration.previous_ksp;
              CPPUNIT_ASSERT_MESSAGE("the SNES does not hold the hook's replacement KSP",
                                     passes_everywhere(previous_ksp != nullptr &&
                                                       ksp != previous_ksp));

              // The first solve's KSP still holds the PC the options step retyped to
              // ILU. A solver that read the KSP or the PC before the hooks would have
              // forced PCNONE on it; one that reads them after the hooks never touches it.
              PC previous_pc;
              LibmeshPetscCallA(TestCommWorld->get(), KSPGetPC(previous_ksp, &previous_pc));
              CPPUNIT_ASSERT_MESSAGE("the previous KSP holds the current PC",
                                     passes_everywhere(previous_pc != pc));
              PetscBool is_none = PETSC_FALSE;
              LibmeshPetscCallA(TestCommWorld->get(),
                                PetscObjectTypeCompare((PetscObject)previous_pc, PCNONE, &is_none));
              CPPUNIT_ASSERT_MESSAGE(
                "the solver forced PCNONE on a PC that was no longer installed",
                passes_everywhere(is_none == PETSC_FALSE));
            }
        }

      CPPUNIT_ASSERT_MESSAGE("the configuration hook did not run on both solves",
                             passes_everywhere(configuration.n_configures == 2));
    }
  }

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Case 5: the matrix-free mode is turned off between two solves of one retained,
  // already set-up SNES. Nothing can rebuild its operators for the new mode, so the
  // second solve must say so instead of silently keeping the operator of the first.
  // The options are set to "false" rather than cleared: PETSc writes a flag only when
  // its option is present, so a cleared option leaves the SNES matrix-free.
  void testReuseRejectsModeChange()
  {
    LOG_UNIT_TEST;

    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf"});
    set_matrix_free_options(options, true);

    {
      CubicResidual cubic_residual;
      CubicJacobian cubic_jacobian;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &cubic_residual;
      system.nonlinear_solver->jacobian_object = &cubic_jacobian;

      system.solve();

      options.set("-snes_mf_operator", "false");
      options.set("-snes_mf", "false");

      // The error must be the solver's own, which names the way out.
      bool rejected = false;
      bool names_the_way_out = false;
      try
        {
          system.solve();
        }
      catch (const libMesh::LogicError & e)
        {
          rejected = true;
          names_the_way_out =
            std::string(e.what()).find("force_new_preconditioner()") != std::string::npos;
        }
      CPPUNIT_ASSERT_MESSAGE("a matrix-free mode changed on a retained SNES was not rejected",
                             passes_everywhere(rejected));
      CPPUNIT_ASSERT_MESSAGE("the error does not name force_new_preconditioner()",
                             passes_everywhere(names_the_way_out));
    }
  }
#endif

  // Case 6: -snes_fd beside -snes_mf, with a Jacobian object registered. The options step
  // installs PETSc's finite-difference Jacobian callback on every solve, and on solve 1 the
  // matrix-free setup that follows replaces it with the matrix-free one. On a retained SNES
  // that setup does not run, so the solver has to put the matrix-free callback back after
  // the options step: put back before it, the finite-difference callback is what solve 2
  // runs, and it tries to assemble a Jacobian into the matrix-free preconditioning matrix.
  void testReuseKeepsMatrixFreeUnderFiniteDifferenceOption()
  {
    LOG_UNIT_TEST;

    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf", "-snes_fd"});
    set_matrix_free_options(options, false);
    options.set("-snes_fd", "true");

    {
      CubicResidual cubic_residual;
      CubicJacobian cubic_jacobian;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &cubic_residual;
      system.nonlinear_solver->jacobian_object = &cubic_jacobian;

      for (unsigned int solve_number = 1; solve_number <= 2; ++solve_number)
        {
          system.solution->zero(); // away from the root: see case 1
          system.solution->close();
          system.update();
          system.solve();
          CPPUNIT_ASSERT_MESSAGE("the solve did no Newton iteration",
                                 passes_everywhere(system.n_nonlinear_iterations() > 0));

          SNESConvergedReason reason;
          LibmeshPetscCallA(TestCommWorld->get(),
                            SNESGetConvergedReason(petsc_solver(system).snes(), &reason));
          CPPUNIT_ASSERT_MESSAGE("the matrix-free solve with -snes_fd did not converge",
                                 passes_everywhere(reason > 0));
          CPPUNIT_ASSERT_MESSAGE("the solution is not the root of x^3 + x = 1",
                                 solution_is_cubic_root(system));

          Mat Amat, Pmat;
          get_matrices(petsc_solver(system), Amat, Pmat);
          PetscBool is_mffd = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)Amat, MATMFFD, &is_mffd));
          CPPUNIT_ASSERT_MESSAGE("SNES operator is not matrix-free",
                                 passes_everywhere(is_mffd == PETSC_TRUE));
          CPPUNIT_ASSERT_MESSAGE("-snes_mf must use one matrix for both",
                                 passes_everywhere(Amat == Pmat));
        }
    }
  }

#if defined(LIBMESH_ENABLE_EXCEPTIONS) && !defined(PETSC_USE_DEBUG)
  // Case 7: a left nonlinear preconditioner, with options that stay the same between the
  // two solves of one retained SNES. Under a left nonlinear preconditioner PETSc's
  // SNESSetUp() rewrites the flags to mf = true, mf_operator = false, and the options step
  // of the next solve reads -snes_mf_operator back, so a solver that compared the two
  // would reject a mode change that never happened. The system is residual only, so the
  // operator is the matrix-free one that setup installs, and both solves start at the root
  // of R(x) = x: they converge at the initial iterate, before the nonlinear preconditioner
  // or a Jacobian is needed. Against a debug PETSc that setup stops earlier, in the check
  // of DMCreateMatrix() that the DM attached itself to the matrix it returned, which
  // libMesh's DM does not do; the case is left out there.
  void testReuseWithLeftNonlinearPreconditioner()
  {
    LOG_UNIT_TEST;

    PetscOptionsScope options({"-snes_mf_operator", "-snes_mf", "-snes_type", "-snes_npc_side",
                               "-npc_snes_type", "-snes_function_type"});
    set_matrix_free_options(options, true);
    options.set("-snes_type", "newtonls");
    options.set("-snes_npc_side", "left");
    options.set("-npc_snes_type", "newtonls");
    options.set("-snes_function_type", "unpreconditioned");

    {
      IdentityResidual identity_residual;
      Mesh mesh(*TestCommWorld);
      EquationSystems es(mesh);
      NonlinearImplicitSystem & system = build(es, mesh);

      system.nonlinear_solver->residual_object = &identity_residual;

      for (unsigned int solve_number = 1; solve_number <= 2; ++solve_number)
        {
          system.solution->zero(); // the root of R(x) = x
          system.solution->close();
          system.update();

          // Solve 2 is the one that regresses: the solver's mode check rejects it, with the
          // error that names force_new_preconditioner(). Any other error propagates.
          bool rejected = false;
          try
            {
              system.solve();
            }
          catch (const libMesh::LogicError & e)
            {
              if (std::string(e.what()).find("force_new_preconditioner()") == std::string::npos)
                throw;
              rejected = true;
            }
          CPPUNIT_ASSERT_MESSAGE("an unchanged matrix-free mode was rejected",
                                 passes_everywhere(!rejected));

          SNESConvergedReason reason;
          LibmeshPetscCallA(TestCommWorld->get(),
                            SNESGetConvergedReason(petsc_solver(system).snes(), &reason));
          CPPUNIT_ASSERT_MESSAGE("the solve did not converge",
                                 passes_everywhere(reason > 0));
          CPPUNIT_ASSERT_MESSAGE("the solve did not stop at its initial iterate",
                                 passes_everywhere(system.n_nonlinear_iterations() == 0));

          // Without a matrix-free operator solve 2 would not reach the mode check at all.
          Mat Amat, Pmat;
          get_matrices(petsc_solver(system), Amat, Pmat);
          PetscBool is_mffd = PETSC_FALSE;
          LibmeshPetscCallA(TestCommWorld->get(),
                            PetscObjectTypeCompare((PetscObject)Amat, MATMFFD, &is_mffd));
          CPPUNIT_ASSERT_MESSAGE("SNES operator is not matrix-free",
                                 passes_everywhere(is_mffd == PETSC_TRUE));
        }
    }
  }
#endif

  void testReuseKeepsMatrixFreeOperator()      { LOG_UNIT_TEST; runTwoSolves(true); }
  void testReuseKeepsFullyMatrixFreeOperator() { LOG_UNIT_TEST; runTwoSolves(false); }
};

CPPUNIT_TEST_SUITE_REGISTRATION( NonlinearImplicitSystemReuseTest );

#endif // LIBMESH_HAVE_PETSC && !PETSC_VERSION_LESS_THAN(3,8,0)
