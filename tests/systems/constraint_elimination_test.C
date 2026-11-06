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
#include <limits>
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
  /**
   * Builds the finite element objects for the system's single variable
   */
  void init(const NonlinearImplicitSystem & sys)
  {
    if (_fe)
      return;
    const FEType fe_type = sys.get_dof_map().variable_type(0);
    _fe = FEBase::build(2, fe_type);
    _qrule = std::make_unique<QGauss>(2, fe_type.default_quadrature_order());
    _fe->attach_quadrature_rule(_qrule.get());
    _fe->get_JxW();
    _fe->get_phi();
    _fe->get_dphi();
  }

  /**
   * Computes the unconstrained element residual \p Re and Jacobian \p Ke on \p elem, whose degree
   * of freedom values are \p ue
   */
  void element(const Elem & elem, const DenseVector<Number> & ue, DenseVector<Number> & Re,
               DenseMatrix<Number> & Ke)
  {
    _fe->reinit(&elem);
    const auto & JxW = _fe->get_JxW();
    const auto & phi = _fe->get_phi();
    const auto & dphi = _fe->get_dphi();
    const auto n = ue.size();
    Re.resize(n);
    Ke.resize(n, n);
    for (const auto qp : index_range(JxW))
      {
        Number u = 0;
        Gradient grad_u;
        for (const auto j : make_range(n))
          {
            u += phi[j][qp] * ue(j);
            grad_u.add_scaled(dphi[j][qp], ue(j));
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
  }

  virtual void residual_and_jacobian(const NumericVector<Number> & soln,
                                     NumericVector<Number> * residual,
                                     SparseMatrix<Number> * jacobian,
                                     NonlinearImplicitSystem & sys) override
  {
    const DofMap & dof_map = sys.get_dof_map();
    this->init(sys);

    std::vector<dof_id_type> dofs;
    DenseVector<Number> ue, Re;
    DenseMatrix<Number> Ke;
    for (const Elem * elem : sys.get_mesh().active_local_element_ptr_range())
      {
        dof_map.dof_indices(elem, dofs);
        ue.resize(dofs.size());
        for (const auto j : index_range(dofs))
          ue(j) = soln(dofs[j]);
        this->element(*elem, ue, Re, Ke);
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

private:
  std::unique_ptr<FEBase> _fe;
  std::unique_ptr<QBase> _qrule;
};

#ifdef LIBMESH_HAVE_EIGEN
/**
 * Nonlinear elimination of the element interior degrees of freedom of \p NonlinearDiffusion: on
 * each element, Newton's method solves the interior rows of the element residual for the interior
 * values, given the element boundary values. Interior degrees of freedom belong to one element
 * and are never constrained, so each local problem is independent of the others
 */
class InteriorElimination : public NonlinearImplicitSystem::ComputeEliminatedSolution
{
public:
  explicit InteriorElimination(NonlinearDiffusion & assembly) : _assembly(assembly) {}

  virtual bool eliminate(NonlinearImplicitSystem & sys) override
  {
    _assembly.init(sys);
    const DofMap & dof_map = sys.get_dof_map();
    const auto & sc_map = dof_map.get_static_condensation();
    const NumericVector<Number> & soln = *sys.current_local_solution;

    std::vector<dof_id_type> dofs, interior;
    std::vector<unsigned int> interior_local;
    DenseVector<Number> ue, Re, R1, du;
    DenseMatrix<Number> Ke, K11;
    for (const Elem * elem : sys.get_mesh().active_local_element_ptr_range())
      {
        dof_map.dof_indices(elem, dofs);
        ue.resize(dofs.size());
        interior.clear();
        interior_local.clear();
        for (const auto j : index_range(dofs))
          {
            ue(j) = soln(dofs[j]);
            // Condensed degrees of freedom are those the parent space retains and static
            // condensation does not. Constrained degrees of freedom are in neither space
            if (sc_map.parent_dof(dofs[j]) != DofObject::invalid_id &&
                sc_map.find_reduced_dof(dofs[j]) == DofObject::invalid_id)
              {
                interior.push_back(dofs[j]);
                interior_local.push_back(j);
              }
          }
        const auto n1 = interior_local.size();
        if (!n1)
          continue;

        // Newton's method on the interior rows. It converges quadratically, so iterating to
        // roundoff costs only a few iterations more than a loose tolerance and keeps the reduced
        // residual a smooth function of the boundary values
        bool converged = false;
        Real initial_norm = -1;
        for (unsigned int it = 0; it < max_its && !converged; ++it)
          {
            _assembly.element(*elem, ue, Re, Ke);
            R1.resize(n1);
            K11.resize(n1, n1);
            for (const auto a : make_range(n1))
              {
                R1(a) = Re(interior_local[a]);
                for (const auto b : make_range(n1))
                  K11(a, b) = Ke(interior_local[a], interior_local[b]);
              }
            const Real norm = R1.l2_norm();
            if (initial_norm < 0)
              initial_norm = norm;
            if (norm <= tolerance * initial_norm || norm < absolute_tolerance)
              {
                converged = true;
                break;
              }
            K11.lu_solve(R1, du);
            for (const auto a : make_range(n1))
              ue(interior_local[a]) -= du(a);
          }
        if (!converged)
          return false;

        for (const auto a : make_range(n1))
          sys.solution->set(interior[a], ue(interior_local[a]));
        ++n_solves;
      }
    return true;
  }

  /// Relative tolerance on the interior residual, near machine precision
  Real tolerance = 1e2 * std::numeric_limits<Real>::epsilon();
  /// Absolute tolerance on the interior residual, for interior values that are already converged
  Real absolute_tolerance = 1e-14;
  unsigned int max_its = 20;
  /// The number of element solves performed, to check that elimination ran
  unsigned int n_solves = 0;

private:
  NonlinearDiffusion & _assembly;
};
#endif

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
  CPPUNIT_TEST(testNonlinearElimination);
  CPPUNIT_TEST(testNonlinearEliminationEnforcedConstraints);
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
    // A distributed mesh keeps periodic neighbors only if it knows the periodic boundaries before
    // it deletes remote elements
    mesh.allow_remote_element_removal(false);
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

    // Ghost the point neighbors of periodic images, which periodic constraints need on a
    // distributed mesh. Both systems have the same periodic boundaries.
    mesh.default_ghosting().set_periodic_boundaries(
        full_sys.get_dof_map().get_periodic_boundaries());

    // Refine one element to create hanging nodes
    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

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
    // A distributed mesh keeps periodic neighbors only if it knows the periodic boundaries before
    // it deletes remote elements
    mesh.allow_remote_element_removal(false);
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

    // Ghost the point neighbors of periodic images, which periodic constraints need on a
    // distributed mesh. Both systems have the same periodic boundaries.
    mesh.default_ghosting().set_periodic_boundaries(
        full_sys.get_dof_map().get_periodic_boundaries());
    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

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
   * the reduced solve as the bounds of the unconstrained degrees of freedom. Elimination bounds
   * only the unconstrained degrees of freedom, so a hanging node value C y can exceed the bound,
   * and the two solves agree when the bound is inactive at constrained degrees of freedom
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

    // Refine at the Dirichlet side, away from where the unbounded solution peaks, so that the
    // bound is inactive at the hanging nodes. Element ids depend on the partitioning of a
    // distributed mesh, so the element is chosen by position
    for (Elem * elem : mesh.active_element_ptr_range())
      if (const Point c = elem->vertex_average(); c(0) < 0.25 && c(1) > 0.5 && c(1) < 0.75)
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
    // A distributed mesh keeps periodic neighbors only if it knows the periodic boundaries before
    // it deletes remote elements
    mesh.allow_remote_element_removal(false);
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

    // Ghost the point neighbors of periodic images, which periodic constraints need on a
    // distributed mesh
    mesh.default_ghosting().set_periodic_boundaries(sys.get_dof_map().get_periodic_boundaries());
    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

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
    // A distributed mesh keeps periodic neighbors only if it knows the periodic boundaries before
    // it deletes remote elements
    mesh.allow_remote_element_removal(false);
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

    // Ghost the point neighbors of periodic images, which periodic constraints need on a
    // distributed mesh. Both systems have the same periodic boundaries.
    mesh.default_ghosting().set_periodic_boundaries(
        plain_sys.get_dof_map().get_periodic_boundaries());
    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

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
   * Nonlinear elimination of the element interior degrees of freedom: SNES solves for the element
   * boundary degrees of freedom only, and each function evaluation first solves the interior
   * equations element by element. The solution must match a plain Newton solve. The Jacobian of
   * the reduced residual is the static condensation Schur complement at the eliminated state, so
   * it must match a finite difference Jacobian of the reduced residual
   */
  void testNonlinearElimination()
  {
    LOG_UNIT_TEST;
    nonlinearElimination(/*eliminate_constraints=*/true);
  }

  /**
   * Nonlinear elimination of a system that keeps its constrained degrees of freedom. The
   * constraint rows of the reduced residual and the Schur complement are enforced instead
   */
  void testNonlinearEliminationEnforcedConstraints()
  {
    LOG_UNIT_TEST;
    nonlinearElimination(/*eliminate_constraints=*/false);
  }

  void nonlinearElimination(const bool eliminate_constraints)
  {
#if !defined(LIBMESH_PETSC_HAVE_MUMPS) && !defined(LIBMESH_PETSC_HAVE_SUPERLU_DIST)
    // The solves use LU, which needs a parallel direct solver when the matrices are distributed
    if (TestCommWorld->size() > 1)
      return;
#endif

    Mesh mesh(*TestCommWorld);
    mesh.allow_remote_element_removal(false);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & plain_sys = es.add_system<NonlinearImplicitSystem>("plain");
    auto & ne_sys = es.add_system<NonlinearImplicitSystem>("ne");
    if (eliminate_constraints)
      {
        plain_sys.create_constraint_elimination();
        ne_sys.create_constraint_elimination();
      }
    ne_sys.create_static_condensation();

    NonlinearDiffusion plain_assembly, ne_assembly;
    InteriorElimination elimination(ne_assembly);
    for (auto [sys, assembly] : {std::make_pair(&plain_sys, &plain_assembly),
                                 std::make_pair(&ne_sys, &ne_assembly)})
      {
        // Third order hierarchic shape functions have interior dofs to eliminate
        const auto u = sys->add_variable("u", THIRD, HIERARCHIC);
        PeriodicBoundary periodic(RealVectorValue(0., 1.));
        periodic.myboundary = 0;
        periodic.pairedboundary = 2;
        sys->get_dof_map().add_periodic_boundary(periodic);
        sys->get_dof_map().add_dirichlet_boundary(
            DirichletBoundary({3}, {u}, ConstFunction<Number>(2.)));
        sys->nonlinear_solver->residual_and_jacobian_object = assembly;
        sys->nonlinear_solver->relative_residual_tolerance = TOLERANCE * TOLERANCE;
      }
    ne_sys.nonlinear_solver->elimination_object = &elimination;

    mesh.default_ghosting().set_periodic_boundaries(
        plain_sys.get_dof_map().get_periodic_boundaries());
    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

    if (Elem * elem = mesh.query_elem_ptr(5); elem && elem->active())
      elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();
    es.init();

    // Nonlinear elimination solves the reduced system with an ordinary linear solver
    const std::vector<std::pair<const char *, const char *>> options = {
        {"-ne_ksp_type", "preonly"}, {"-ne_pc_type", "lu"}, {"-plain_pc_type", "lu"}};
    for (const auto & [name, value] : options)
      LibmeshPetscCall2(mesh.comm(), PetscOptionsSetValue(nullptr, name, value));
    // Compares the Schur complement with a finite difference Jacobian of the reduced residual at
    // the first Newton iterate, after the interior equations have been solved once
    struct JacobianTest
    {
      static PetscErrorCode monitor(SNES snes, PetscInt its, PetscReal, void * ctx)
      {
        PetscFunctionBeginUser;
        auto & self = *static_cast<JacobianTest *>(ctx);
        if (its == 1)
          {
            // SNESTestJacobian compares the Jacobian SNES holds, so assemble it at this iterate
            Vec x;
            Mat J, P;
            PetscCall(SNESGetSolution(snes, &x));
            PetscCall(SNESGetJacobian(snes, &J, &P, nullptr, nullptr));
            PetscCall(SNESComputeJacobian(snes, x, J, P));
            PetscCall(SNESTestJacobian(snes, &self.jac_norm, &self.diff_norm));
            self.tested = true;
          }
        PetscFunctionReturn(PETSC_SUCCESS);
      }
      PetscReal jac_norm = 0, diff_norm = 0;
      bool tested = false;
    } jacobian_test;

    for (auto * sys : {&plain_sys, &ne_sys})
      {
        sys->prefix_with_name(true);
        sys->solution->zero();
        sys->solution->close();
        sys->update();
        if (sys == &ne_sys)
          {
            auto & solver = cast_ref<PetscNonlinearSolver<Number> &>(*sys->nonlinear_solver);
            LibmeshPetscCall2(mesh.comm(),
                              SNESMonitorSet(solver.snes(sys->prefix().c_str()),
                                             JacobianTest::monitor,
                                             &jacobian_test,
                                             nullptr));
          }
        sys->solve();
        CPPUNIT_ASSERT(sys->nonlinear_solver->converged);
      }
    for (const auto & [name, value] : options)
      {
        libmesh_ignore(value);
        LibmeshPetscCall2(mesh.comm(), PetscOptionsClearValue(nullptr, name));
      }

    CPPUNIT_ASSERT(elimination.n_solves > 0);

    // Finite differences of the reduced residual re-solve the interior equations at each perturbed
    // iterate, so they include the dependence of the interior values on the boundary ones. Forward
    // differences are accurate to about the square root of machine precision
    CPPUNIT_ASSERT(jacobian_test.tested);
    CPPUNIT_ASSERT_LESS(1e-6 * jacobian_test.jac_norm, jacobian_test.diff_norm);

    const auto & dof_map = plain_sys.get_dof_map();
    for (const auto i : make_range(dof_map.first_dof(), dof_map.end_dof()))
      LIBMESH_ASSERT_NUMBERS_EQUAL(
          (*plain_sys.solution)(i), (*ne_sys.solution)(i), std::sqrt(TOLERANCE) * TOLERANCE);
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
