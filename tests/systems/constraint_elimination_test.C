#include <libmesh/constraint_elimination_dof_map.h>
#include <libmesh/const_function.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/nonlinear_solver.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dense_vector.h>
#include <tuple>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/petsc_nonlinear_solver.h>
#include <libmesh/petscdmlibmesh.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/reduced_sparse_matrix.h>
#include <libmesh/static_condensation_dof_map.h>
#include <libmesh/static_condensation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

// The tests use PETSc solvers and a mesh with hanging nodes, periodic and Dirichlet constraints
#if defined(LIBMESH_HAVE_PETSC) && defined(LIBMESH_ENABLE_AMR) &&                                  \
    defined(LIBMESH_ENABLE_PERIODIC) && defined(LIBMESH_ENABLE_DIRICHLET) && LIBMESH_DIM > 1

// Defined in petscdmlibmeshimpl.C without a public declaration
PetscErrorCode DMlibMeshCreateFieldDecompositionDM(DM dm, PetscInt dnumber, PetscInt * dsizes, char *** dvarlists, DM * ddm);

namespace
{
/**
 * Assembles a non-symmetric operator, a Laplacian plus an advection term, so that the test can
 * tell a constrained row of C^T K C from a constrained column
 */
/**
 * Residual and Jacobian of -div((1 + u^2) grad u) + u_x = 1, a nonlinear problem whose Newton
 * iterates pass through states that violate no constraint only if the solver forms C y + h
 */
class NonlinearDiffusion : public NonlinearImplicitSystem::ComputeResidualandJacobian
{
public:
  virtual void residual_and_jacobian(const NumericVector<Number> & soln,
                                     NumericVector<Number> * residual,
                                     SparseMatrix<Number> * jacobian,
                                     NonlinearImplicitSystem & sys) override
  {
    const DofMap & dof_map = sys.get_dof_map();
    const FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe(FEBase::build(2, fe_type));
    QGauss qrule(2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    const auto & JxW = fe->get_JxW();
    const auto & phi = fe->get_phi();
    const auto & dphi = fe->get_dphi();

    std::vector<dof_id_type> dofs;
    DenseVector<Number> Re;
    DenseMatrix<Number> Ke;
    for (const Elem * elem : sys.get_mesh().active_local_element_ptr_range())
      {
        dof_map.dof_indices(elem, dofs);
        fe->reinit(elem);
        const auto n = dofs.size();
        Re.resize(n);
        Ke.resize(n, n);
        for (const auto qp : index_range(JxW))
          {
            Number u = 0;
            Gradient grad_u;
            for (const auto j : make_range(n))
              {
                u += phi[j][qp] * soln(dofs[j]);
                grad_u.add_scaled(dphi[j][qp], soln(dofs[j]));
              }
            const Number k = 1. + u * u;
            for (const auto i : make_range(n))
              {
                Re(i) += JxW[qp] * (k * grad_u * dphi[i][qp] + (grad_u(0) - 1.) * phi[i][qp]);
                for (const auto j : make_range(n))
                  Ke(i, j) += JxW[qp] * (k * dphi[j][qp] * dphi[i][qp] +
                                         2. * u * phi[j][qp] * (grad_u * dphi[i][qp]) +
                                         dphi[j][qp](0) * phi[i][qp]);
              }
          }
        dof_map.constrain_element_matrix_and_vector(Ke, Re, dofs, false);
        if (residual)
          residual->add_vector(Re, dofs);
        if (jacobian)
          {
#ifdef LIBMESH_HAVE_EIGEN
            // A statically condensed Jacobian accumulates element matrices by element
            if (auto * sc = dynamic_cast<StaticCondensation *>(jacobian))
              sc->set_current_elem(*elem);
#endif
            jacobian->add_matrix(Ke, dofs);
          }
      }
  }
};

/**
 * Caps the solution from above, so a variational inequality solve has active bounds
 */
class UpperBound : public NonlinearImplicitSystem::ComputeBounds
{
public:
  explicit UpperBound(Real cap) : _cap(cap) {}

  virtual void bounds(NumericVector<Number> &,
                      NumericVector<Number> & XU,
                      NonlinearImplicitSystem &) override
  {
    for (const auto i : make_range(XU.first_local_index(), XU.last_local_index()))
      XU.set(i, _cap);
    XU.close();
  }

private:
  const Real _cap;
};

/**
 * The residual and Jacobian of the Laplacian with a zero-mean source, -div(grad u) = x - 1/2, whose
 * fully periodic solutions differ by constants
 */
class PeriodicPoisson : public NonlinearImplicitSystem::ComputeResidualandJacobian
{
public:
  virtual void residual_and_jacobian(const NumericVector<Number> & soln,
                                     NumericVector<Number> * residual,
                                     SparseMatrix<Number> * jacobian,
                                     NonlinearImplicitSystem & sys) override
  {
    const DofMap & dof_map = sys.get_dof_map();
    const FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe(FEBase::build(2, fe_type));
    QGauss qrule(2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    const auto & JxW = fe->get_JxW();
    const auto & phi = fe->get_phi();
    const auto & dphi = fe->get_dphi();
    const auto & xyz = fe->get_xyz();

    std::vector<dof_id_type> dofs;
    DenseVector<Number> Re;
    DenseMatrix<Number> Ke;
    for (const Elem * elem : sys.get_mesh().active_local_element_ptr_range())
      {
        dof_map.dof_indices(elem, dofs);
        fe->reinit(elem);
        const auto n = dofs.size();
        Re.resize(n);
        Ke.resize(n, n);
        for (const auto qp : index_range(JxW))
          for (const auto i : make_range(n))
            {
              Re(i) -= JxW[qp] * (xyz[qp](0) - 0.5) * phi[i][qp];
              for (const auto j : make_range(n))
                {
                  Ke(i, j) += JxW[qp] * dphi[j][qp] * dphi[i][qp];
                  Re(i) += JxW[qp] * dphi[j][qp] * dphi[i][qp] * soln(dofs[j]);
                }
            }
        dof_map.constrain_element_matrix_and_vector(Ke, Re, dofs, false);
        if (residual)
          residual->add_vector(Re, dofs);
        if (jacobian)
          jacobian->add_matrix(Ke, dofs);
      }
  }
};

/**
 * The constants, the null space of the fully periodic Laplacian
 */
class Constants : public NonlinearImplicitSystem::ComputeVectorSubspace
{
public:
  virtual void operator()(std::vector<NumericVector<Number> *> & sp,
                          NonlinearImplicitSystem & sys) override
  {
    _mode = sys.solution->zero_clone();
    _mode->add(1.);
    _mode->close();
    sp.push_back(_mode.get());
  }

private:
  std::unique_ptr<NumericVector<Number>> _mode;
};

void assemble_operator(const System & sys, SparseMatrix<Number> & mat)
{
  const MeshBase & mesh = sys.get_mesh();
  const DofMap & dof_map = sys.get_dof_map();
  const FEType fe_type = dof_map.variable_type(0);

  std::unique_ptr<FEBase> fe(FEBase::build(2, fe_type));
  QGauss qrule(2, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const auto & JxW = fe->get_JxW();
  const auto & phi = fe->get_phi();
  const auto & dphi = fe->get_dphi();

  std::vector<dof_id_type> dof_indices;
  DenseMatrix<Number> Ke;

  for (const Elem * elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices(elem, dof_indices);
      fe->reinit(elem);
      Ke.resize(dof_indices.size(), dof_indices.size());
      for (const auto qp : index_range(JxW))
        for (const auto i : index_range(dof_indices))
          for (const auto j : index_range(dof_indices))
            Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp] + phi[i][qp] * dphi[j][qp](0));

      // Homogeneous condensation: the element matrix becomes C^T K C on unconstrained rows and
      // columns
      dof_map.constrain_element_matrix(Ke, dof_indices, false);
      mat.add_matrix(Ke, dof_indices);
    }
  mat.close();
}
}

class ConstraintEliminationTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(ConstraintEliminationTest);

  CPPUNIT_TEST(testReducedOperator);
  CPPUNIT_TEST(testReducedSolve);
  CPPUNIT_TEST(testReducedFieldDecomposition);
  CPPUNIT_TEST(testReducedBounds);
  CPPUNIT_TEST(testReducedNullSpace);
#ifdef LIBMESH_HAVE_EIGEN
  CPPUNIT_TEST(testCondensedReducedSolve);
#endif

  CPPUNIT_TEST_SUITE_END();

private:
  /**
   * On a mesh with hanging nodes, a periodic pair of sides and a Dirichlet side, the reduced
   * matrix assembled with full indices must equal the unconstrained rows and columns of the full
   * constrained matrix, and every constrained degree of freedom must be eliminated
   */
  void testReducedOperator()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    // build_square's boundary ids: 0 bottom, 1 right, 2 top, 3 left
    EquationSystems es(mesh);
    auto & full_sys = es.add_system<LinearImplicitSystem>("full");
    auto & reduced_sys = es.add_system<LinearImplicitSystem>("reduced_sys");
    reduced_sys.create_constraint_elimination();

    for (System * sys : {static_cast<System *>(&full_sys), static_cast<System *>(&reduced_sys)})
      {
        const auto u = sys->add_variable("u", SECOND, LAGRANGE);
        PeriodicBoundary periodic(RealVectorValue(0., 1.));
        periodic.myboundary = 0;
        periodic.pairedboundary = 2;
        sys->get_dof_map().add_periodic_boundary(periodic);
        sys->get_dof_map().add_dirichlet_boundary(
            DirichletBoundary({3}, {u}, ConstFunction<Number>(1.)));
      }

    // Refine one element to create hanging nodes
    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();

    es.init();

    const DofMap & dof_map = full_sys.get_dof_map();
    const auto & reduced_map = reduced_sys.get_dof_map().get_constraint_elimination();
    CPPUNIT_ASSERT(dof_map.n_constrained_dofs() > 0);
    CPPUNIT_ASSERT_EQUAL(dof_map.n_dofs() - dof_map.n_constrained_dofs(), reduced_map.n_dofs());

    assemble_operator(full_sys, *full_sys.matrix);
    assemble_operator(reduced_sys, *reduced_sys.matrix);

    // Both systems number their full dofs alike, so the reduced matrix's rows and columns are
    // the unconstrained rows and columns of the full matrix
    const std::vector<numeric_index_type> retained(reduced_map.local_retained_dofs().begin(),
                                                   reduced_map.local_retained_dofs().end());
    auto expected = SparseMatrix<Number>::build(full_sys.comm());
    full_sys.matrix->create_submatrix(*expected, retained, retained);

    const auto & reduced = cast_ref<ReducedSparseMatrix<Number> &>(*reduced_sys.matrix);
    std::vector<numeric_index_type> cols;
    std::vector<Number> vals;
    for (const auto i : index_range(retained))
      {
        expected->get_row(expected->row_start() + i, cols, vals);
        for (const auto k : index_range(cols))
          LIBMESH_ASSERT_NUMBERS_EQUAL(
              vals[k], reduced(reduced_map.first_dof() + i, cols[k]), TOLERANCE * TOLERANCE);
      }

    // The reduced operator has no other entries
    LIBMESH_ASSERT_FP_EQUAL(expected->l1_norm(), reduced.l1_norm(), TOLERANCE * TOLERANCE);
  }

  /**
   * Solves a nonlinear problem on the same constrained mesh with and without eliminating the
   * constrained degrees of freedom. Both must converge to the same solution, and eliminating must
   * keep every iterate on the constraints, including a nonzero Dirichlet value
   */
  void testReducedSolve()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & full_sys = es.add_system<NonlinearImplicitSystem>("full");
    auto & reduced_sys = es.add_system<NonlinearImplicitSystem>("reduced_sys");
    reduced_sys.create_constraint_elimination();

    NonlinearDiffusion assembly;
    for (auto * sys : {&full_sys, &reduced_sys})
      {
        const auto u = sys->add_variable("u", SECOND, LAGRANGE);
        PeriodicBoundary periodic(RealVectorValue(0., 1.));
        periodic.myboundary = 0;
        periodic.pairedboundary = 2;
        sys->get_dof_map().add_periodic_boundary(periodic);
        sys->get_dof_map().add_dirichlet_boundary(
            DirichletBoundary({3}, {u}, ConstFunction<Number>(2.)));
        sys->nonlinear_solver->residual_and_jacobian_object = &assembly;
        sys->nonlinear_solver->relative_residual_tolerance = TOLERANCE * TOLERANCE;
      }

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();

    es.init();

    for (auto * sys : {&full_sys, &reduced_sys})
      {
        // Start away from the Dirichlet value, so constrained entries must be imposed
        sys->solution->zero();
        sys->solution->close();
        sys->update();
        sys->solve();
        CPPUNIT_ASSERT(sys->nonlinear_solver->converged);
      }

    const auto & dof_map = full_sys.get_dof_map();
    for (const auto i : make_range(dof_map.first_dof(), dof_map.end_dof()))
      LIBMESH_ASSERT_NUMBERS_EQUAL(
          (*full_sys.solution)(i), (*reduced_sys.solution)(i), std::sqrt(TOLERANCE) * TOLERANCE);
  }

  /**
   * A variational inequality solve with an active upper bound must give the same solution with
   * and without eliminating constraints, so the bounds a user computes on full vectors must reach
   * the reduced solve as the bounds of the unconstrained degrees of freedom
   */
  void testReducedBounds()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & full_sys = es.add_system<NonlinearImplicitSystem>("full");
    auto & reduced_sys = es.add_system<NonlinearImplicitSystem>("reduced_sys");
    reduced_sys.create_constraint_elimination();

    // The unbounded solution peaks near 1.2, above this cap
    const Real cap_value = 1.1;
    NonlinearDiffusion assembly;
    UpperBound cap(cap_value);
    for (auto * sys : {&full_sys, &reduced_sys})
      {
        const auto u = sys->add_variable("u", SECOND, LAGRANGE);
        sys->get_dof_map().add_dirichlet_boundary(
            DirichletBoundary({3}, {u}, ConstFunction<Number>(1.)));
        sys->nonlinear_solver->residual_and_jacobian_object = &assembly;
        sys->nonlinear_solver->bounds_object = &cap;
        sys->nonlinear_solver->relative_residual_tolerance = TOLERANCE * TOLERANCE;
      }

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();
    es.init();

    LibmeshPetscCall2(mesh.comm(), PetscOptionsSetValue(nullptr, "-snes_type", "vinewtonrsls"));
    for (auto * sys : {&full_sys, &reduced_sys})
      {
        sys->solution->zero();
        sys->solution->close();
        sys->update();
        sys->solve();
        CPPUNIT_ASSERT(sys->nonlinear_solver->converged);
      }
    LibmeshPetscCall2(mesh.comm(), PetscOptionsClearValue(nullptr, "-snes_type"));

    const auto & dof_map = full_sys.get_dof_map();
    Real max_u = 0;
    for (const auto i : make_range(dof_map.first_dof(), dof_map.end_dof()))
      {
        LIBMESH_ASSERT_NUMBERS_EQUAL(
            (*full_sys.solution)(i), (*reduced_sys.solution)(i), std::sqrt(TOLERANCE) * TOLERANCE);
        max_u = std::max(max_u, std::abs((*reduced_sys.solution)(i)));
      }
    mesh.comm().max(max_u);
    // The bound is active
    LIBMESH_ASSERT_FP_EQUAL(cap_value, max_u, TOLERANCE);
  }

  /**
   * A fully periodic Laplacian on a mesh with hanging nodes is singular, with the constants as its
   * null space. A user computes that null space on full vectors; the reduced Jacobian must carry
   * its restriction, which the reduced Jacobian annihilates, and the solve must converge
   */
  void testReducedNullSpace()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & sys = es.add_system<NonlinearImplicitSystem>("sys");
    sys.create_constraint_elimination();
    sys.add_variable("u", SECOND, LAGRANGE);
    for (const auto & [mine, paired, offset] :
         {std::make_tuple(0, 2, RealVectorValue(0., 1.)),
          std::make_tuple(3, 1, RealVectorValue(1., 0.))})
      {
        PeriodicBoundary periodic(offset);
        periodic.myboundary = mine;
        periodic.pairedboundary = paired;
        sys.get_dof_map().add_periodic_boundary(periodic);
      }

    PeriodicPoisson assembly;
    Constants constants;
    sys.nonlinear_solver->residual_and_jacobian_object = &assembly;
    sys.nonlinear_solver->nullspace_object = &constants;
    sys.nonlinear_solver->relative_residual_tolerance = TOLERANCE * TOLERANCE;

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();
    es.init();

    CPPUNIT_ASSERT(sys.get_dof_map().n_constrained_dofs() > 0);

    sys.solution->zero();
    sys.solution->close();
    sys.update();
    sys.solve();
    CPPUNIT_ASSERT(sys.nonlinear_solver->converged);

    // The solver attaches the null space to the system matrix, which is the reduced Jacobian
    const auto & reduced = cast_ref<ReducedSparseMatrix<Number> &>(*sys.matrix);
    MatNullSpace nullspace;
    LibmeshPetscCall2(mesh.comm(), MatGetNullSpace(reduced.mat(), &nullspace));
    CPPUNIT_ASSERT(nullspace);
    PetscBool is_null;
    LibmeshPetscCall2(mesh.comm(), MatNullSpaceTest(nullspace, reduced.mat(), &is_null));
    CPPUNIT_ASSERT(is_null);
  }

  /**
   * Static condensation nested in constraint elimination: the outer Krylov solve acts on the
   * unconstrained degrees of freedom and the condensed solve on their element boundary subset.
   * Condensing exactly must reproduce the eliminated solve without condensation
   */
  void testCondensedReducedSolve()
  {
    LOG_UNIT_TEST;

#if !defined(LIBMESH_PETSC_HAVE_MUMPS) && !defined(LIBMESH_PETSC_HAVE_SUPERLU_DIST)
    // The solves use LU, which needs a parallel direct solver when the matrices are distributed
    if (TestCommWorld->size() > 1)
      return;
#endif

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & plain_sys = es.add_system<NonlinearImplicitSystem>("plain");
    auto & condensed_sys = es.add_system<NonlinearImplicitSystem>("condensed");
    plain_sys.create_constraint_elimination();
    condensed_sys.create_constraint_elimination();
    condensed_sys.create_static_condensation();

    NonlinearDiffusion assembly;
    for (auto * sys : {&plain_sys, &condensed_sys})
      {
        // Third order hierarchic shape functions have interior dofs to condense
        const auto u = sys->add_variable("u", THIRD, HIERARCHIC);
        PeriodicBoundary periodic(RealVectorValue(0., 1.));
        periodic.myboundary = 0;
        periodic.pairedboundary = 2;
        sys->get_dof_map().add_periodic_boundary(periodic);
        sys->get_dof_map().add_dirichlet_boundary(
            DirichletBoundary({3}, {u}, ConstFunction<Number>(2.)));
        sys->nonlinear_solver->residual_and_jacobian_object = &assembly;
        sys->nonlinear_solver->relative_residual_tolerance = TOLERANCE * TOLERANCE;
      }

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();
    es.init();

    // The condensed problem is solved exactly, so the outer Krylov solve needs one iteration. The
    // statically condensed operator has no action of its own, which the default line search needs
    LibmeshPetscCall2(mesh.comm(), PetscOptionsSetValue(nullptr, "-condensed_ksp_type", "preonly"));
    LibmeshPetscCall2(mesh.comm(),
                      PetscOptionsSetValue(nullptr, "-condensed_snes_linesearch_type", "basic"));
    LibmeshPetscCall2(mesh.comm(), PetscOptionsSetValue(nullptr, "-condensed_condensed_pc_type", "lu"));
    LibmeshPetscCall2(mesh.comm(), PetscOptionsSetValue(nullptr, "-plain_pc_type", "lu"));
    for (auto * sys : {&plain_sys, &condensed_sys})
      {
        sys->prefix_with_name(true);
        sys->solution->zero();
        sys->solution->close();
        sys->update();
        sys->solve();
        CPPUNIT_ASSERT(sys->nonlinear_solver->converged);
      }
    for (const char * option : {"-condensed_ksp_type",
                                "-condensed_snes_linesearch_type",
                                "-condensed_condensed_pc_type",
                                "-plain_pc_type"})
      LibmeshPetscCall2(mesh.comm(), PetscOptionsClearValue(nullptr, option));

    const auto & dof_map = plain_sys.get_dof_map();
    CPPUNIT_ASSERT(condensed_sys.get_dof_map().get_static_condensation().n_dofs() <
                   dof_map.get_constraint_elimination().n_dofs());
    for (const auto i : make_range(dof_map.first_dof(), dof_map.end_dof()))
      LIBMESH_ASSERT_NUMBERS_EQUAL(
          (*plain_sys.solution)(i), (*condensed_sys.solution)(i), std::sqrt(TOLERANCE) * TOLERANCE);
  }

  /**
   * A field decomposition of a system that eliminates constraints must partition the reduced space
   * by variable: every local reduced index in exactly one field, each field holding the reduced
   * indices of that variable's unconstrained degrees of freedom
   */
  void testReducedFieldDecomposition()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & sys = es.add_system<NonlinearImplicitSystem>("sys");
    sys.create_constraint_elimination();
    const auto u = sys.add_variable("u", SECOND, LAGRANGE);
    const auto v = sys.add_variable("v", FIRST, LAGRANGE);
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({3}, {u}, ConstFunction<Number>(1.)));
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({1}, {v}, ConstFunction<Number>(0.)));

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();
    es.init();

    auto & solver = cast_ref<PetscNonlinearSolver<Number> &>(*sys.nonlinear_solver);
    DM dm;
    LibmeshPetscCall2(sys.comm(), SNESGetDM(solver.snes(), &dm));
    DM fdm;
    LibmeshPetscCall2(sys.comm(), DMlibMeshCreateFieldDecompositionDM(dm, 0, nullptr, nullptr, &fdm));

    PetscInt n_fields;
    IS * fields;
    LibmeshPetscCall2(sys.comm(),
                      DMCreateFieldDecomposition(fdm, &n_fields, nullptr, &fields, nullptr));
    CPPUNIT_ASSERT_EQUAL(PetscInt(2), n_fields);

    // The expected reduced indices of each variable, from the full numbering
    const DofMap & dof_map = sys.get_dof_map();
    const auto & reduced_map = dof_map.get_constraint_elimination();
    std::vector<std::set<PetscInt>> expected(2);
    std::vector<dof_id_type> dofs;
    for (const unsigned int var : {u, v})
      for (const Elem * elem : mesh.active_local_element_ptr_range())
        {
          dof_map.dof_indices(elem, dofs, var);
          for (const auto dof : dofs)
            if (dof_map.local_index(dof) && !dof_map.is_constrained_dof(dof))
              expected[var].insert(cast_int<PetscInt>(reduced_map.reduced_dof(dof)));
        }

    PetscInt n_found = 0;
    for (const auto f : make_range(n_fields))
      {
        PetscInt n;
        const PetscInt * indices;
        LibmeshPetscCall2(sys.comm(), ISGetLocalSize(fields[f], &n));
        LibmeshPetscCall2(sys.comm(), ISGetIndices(fields[f], &indices));
        const std::set<PetscInt> found(indices, indices + n);
        LibmeshPetscCall2(sys.comm(), ISRestoreIndices(fields[f], &indices));
        CPPUNIT_ASSERT(found == expected[f]);
        n_found += n;
        LibmeshPetscCall2(sys.comm(), ISDestroy(&fields[f]));
      }
    LibmeshPetscCall2(sys.comm(), PetscFree(fields));
    LibmeshPetscCall2(sys.comm(), DMDestroy(&fdm));

    // Together the fields cover the local reduced space once
    CPPUNIT_ASSERT_EQUAL(cast_int<PetscInt>(reduced_map.n_local_dofs()), n_found);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ConstraintEliminationTest);

#endif
