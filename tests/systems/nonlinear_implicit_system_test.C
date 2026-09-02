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

// NonlinearSolver<T>::build() only has a working implementation when built
// against one of these packages; without either, constructing a
// NonlinearImplicitSystem throws unconditionally.
#if defined(LIBMESH_HAVE_PETSC) || defined(LIBMESH_TRILINOS_HAVE_NOX)

#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/nonlinear_solver.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/libmesh_exceptions.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
// Only overrides the single-matrix solve(), like a backend that doesn't
// distinguish the Jacobian operator matrix from the preconditioning matrix,
// relying on NonlinearSolver's own default two-matrix forwarding.
class SingleMatrixNonlinearSolver : public NonlinearSolver<Number>
{
public:
  explicit SingleMatrixNonlinearSolver(NonlinearImplicitSystem & s) : NonlinearSolver<Number>(s) {}

  // Un-hide the base class's two-matrix solve() overload, which we don't
  // override here since it's exactly the default forwarding behavior we're
  // testing.
  using NonlinearSolver<Number>::solve;

  virtual void init(const char * = nullptr) override {}

  virtual std::pair<unsigned int, Real>
  solve(SparseMatrix<Number> & jac_in,
        NumericVector<Number> &,
        NumericVector<Number> &,
        const double,
        const unsigned int) override
  {
    solved_with = &jac_in;
    ++single_matrix_calls;
    return {0, 0.};
  }

  virtual int get_total_linear_iterations() override { return 0; }
  virtual unsigned get_current_nonlinear_iteration_number() const override { return 0; }

  SparseMatrix<Number> * solved_with = nullptr;
  unsigned int single_matrix_calls = 0;
};

// Also records which matrix(es) NonlinearImplicitSystem::solve() actually
// passes down to the NonlinearSolver, without performing any real solve. On
// top of the inherited single-matrix override, this overrides the
// two-matrix solve() overload too, like a backend that distinguishes the
// Jacobian operator matrix from the preconditioning matrix.
class TwoMatrixNonlinearSolver : public SingleMatrixNonlinearSolver
{
public:
  explicit TwoMatrixNonlinearSolver(NonlinearImplicitSystem & s) : SingleMatrixNonlinearSolver(s) {}

  // Un-hide the inherited single-matrix override alongside our own
  // two-matrix override below.
  using SingleMatrixNonlinearSolver::solve;

  virtual std::pair<unsigned int, Real>
  solve(SparseMatrix<Number> & jac_in,
        SparseMatrix<Number> & pre_in,
        NumericVector<Number> &,
        NumericVector<Number> &,
        const double,
        const unsigned int) override
  {
    solved_with = &jac_in;
    solved_with_pre = &pre_in;
    ++two_matrix_calls;
    return {0, 0.};
  }

  SparseMatrix<Number> * solved_with_pre = nullptr;
  unsigned int two_matrix_calls = 0;
};
}

class NonlinearImplicitSystemOperatorMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(NonlinearImplicitSystemOperatorMatrixTest);

  CPPUNIT_TEST(testDefaultUsesSingleMatrixSolve);
  CPPUNIT_TEST(testOperatorMatrixUsesTwoMatrixSolve);
  CPPUNIT_TEST(testTwoMatrixSolveWithoutOverrideThrows);

  CPPUNIT_TEST_SUITE_END();

public:

  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<EquationSystems> es;
  NonlinearImplicitSystem * sys = nullptr;

  void setUp()
  {
    mesh = std::make_unique<Mesh>(*TestCommWorld);
    MeshTools::Generation::build_point(*mesh);
    es = std::make_unique<EquationSystems>(*mesh);
    sys = &es->add_system<NonlinearImplicitSystem>("test");
    es->init();
  }

  void tearDown()
  {
    sys = nullptr;
    es.reset();
    mesh.reset();
  }

  // Without an operator matrix registered, solve() should use the ordinary
  // single-matrix path with the system's own preconditioning matrix.
  void testDefaultUsesSingleMatrixSolve()
  {
    LOG_UNIT_TEST;

    auto solver = std::make_unique<TwoMatrixNonlinearSolver>(*sys);
    auto * solver_ptr = solver.get();
    sys->nonlinear_solver = std::move(solver);

    sys->solve();

    CPPUNIT_ASSERT(solver_ptr->solved_with == sys->matrix);
    CPPUNIT_ASSERT_EQUAL(solver_ptr->single_matrix_calls, 1u);
    CPPUNIT_ASSERT_EQUAL(solver_ptr->two_matrix_calls, 0u);
  }

  // set_operator_matrix() should make solve() use the two-matrix path, with
  // the registered matrix as the Jacobian operator and the system's own
  // matrix as the preconditioning matrix.
  void testOperatorMatrixUsesTwoMatrixSolve()
  {
    LOG_UNIT_TEST;

    auto amat = SparseMatrix<Number>::build(*TestCommWorld);
    sys->set_operator_matrix(amat.get());

    auto solver = std::make_unique<TwoMatrixNonlinearSolver>(*sys);
    auto * solver_ptr = solver.get();
    sys->nonlinear_solver = std::move(solver);

    sys->solve();

    CPPUNIT_ASSERT(solver_ptr->solved_with == amat.get());
    CPPUNIT_ASSERT(solver_ptr->solved_with_pre == sys->matrix);
    CPPUNIT_ASSERT_EQUAL(solver_ptr->two_matrix_calls, 1u);
    CPPUNIT_ASSERT_EQUAL(solver_ptr->single_matrix_calls, 0u);
  }

  // NonlinearSolver's own two-matrix solve() overload has no fallback
  // implementation -- a backend (like this test's mock) that doesn't
  // override it can't honor a distinct Jacobian operator matrix, so it must
  // refuse rather than silently solving with the preconditioning matrix
  // instead of the one the caller explicitly asked to use.
  void testTwoMatrixSolveWithoutOverrideThrows()
  {
    LOG_UNIT_TEST;

    auto amat = SparseMatrix<Number>::build(*TestCommWorld);
    sys->set_operator_matrix(amat.get());

    auto solver = std::make_unique<SingleMatrixNonlinearSolver>(*sys);
    sys->nonlinear_solver = std::move(solver);

    CPPUNIT_ASSERT_THROW(sys->solve(), libMesh::NotImplemented);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(NonlinearImplicitSystemOperatorMatrixTest);

#endif // LIBMESH_HAVE_PETSC || LIBMESH_TRILINOS_HAVE_NOX
