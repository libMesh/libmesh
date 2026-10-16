#include <libmesh/const_function.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/nonlinear_solver.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
/**
 * Residual and Jacobian of -u'' = 0 on [0, 1]. The element contributions
 * are condensed homogeneously, as a NonlinearImplicitSystem user does, so
 * every constrained row, and with it every prescribed value, comes from the
 * solver's own constraint enforcement.
 */
class LaplaceResidualAndJacobian : public NonlinearImplicitSystem::ComputeResidualandJacobian
{
public:
  virtual void residual_and_jacobian(const NumericVector<Number> & soln,
                                     NumericVector<Number> * residual,
                                     SparseMatrix<Number> * jacobian,
                                     NonlinearImplicitSystem & sys) override
  {
    const MeshBase & mesh = sys.get_mesh();
    DofMap & dof_map = sys.get_dof_map();
    const FEType fe_type = dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(1, fe_type));
    QGauss qrule(1, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    const auto & JxW = fe->get_JxW();
    const auto & dphi = fe->get_dphi();

    if (residual)
      residual->zero();
    if (jacobian)
      jacobian->zero();

    std::vector<dof_id_type> dof_indices;
    DenseVector<Number> Re;
    DenseMatrix<Number> Ke;

    for (const Elem * elem : mesh.active_local_element_ptr_range())
      {
        dof_map.dof_indices(elem, dof_indices);
        fe->reinit(elem);

        const auto n = dof_indices.size();
        Re.resize(n);
        Ke.resize(n, n);

        for (auto qp : index_range(JxW))
          for (auto i : make_range(n))
            for (auto j : make_range(n))
              {
                const Number k = JxW[qp] * dphi[i][qp] * dphi[j][qp];
                Ke(i, j) += k;
                Re(i) += k * soln(dof_indices[j]);
              }

        dof_map.constrain_element_matrix_and_vector(Ke, Re, dof_indices, false);

        if (residual)
          residual->add_vector(Re, dof_indices);
        if (jacobian)
          jacobian->add_matrix(Ke, dof_indices);
      }
  }
};
}

class HeterogeneousNonlinearConstraintTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(HeterogeneousNonlinearConstraintTest);

#if defined(LIBMESH_ENABLE_DIRICHLET) && defined(LIBMESH_HAVE_PETSC) && \
    !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  CPPUNIT_TEST(testNonzeroDirichletValues);
#endif

  CPPUNIT_TEST_SUITE_END();

private:
  /**
   * PetscNonlinearSolver enforces a constraint in each residual it hands
   * PETSc. With nonzero DirichletBoundary values on both ends, the exact
   * solution of -u'' = 0 is the line between them, which a converged
   * solve has to reproduce
   */
  void testNonzeroDirichletValues()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line(mesh, 8, 0., 1., EDGE2);

    EquationSystems es(mesh);
    auto & sys = es.add_system<NonlinearImplicitSystem>("sys");
    const auto u = sys.add_variable("u", FIRST, LAGRANGE);

    // build_line's boundary ids
    const boundary_id_type left_id = 0, right_id = 1;
    const Real left_value = 2, right_value = 5;
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({left_id}, {u}, ConstFunction<Number>(left_value)));
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({right_id}, {u}, ConstFunction<Number>(right_value)));

    LaplaceResidualAndJacobian assembly;
    sys.nonlinear_solver->residual_and_jacobian_object = &assembly;

    es.init();

    // Start from an iterate that satisfies neither boundary value, so
    // the constrained rows have to pull it to them
    sys.solution->zero();
    sys.solution->close();
    sys.update();

    sys.solve();

    CPPUNIT_ASSERT(sys.nonlinear_solver->converged);

    for (const Node * node : mesh.local_node_ptr_range())
      {
        const Real x = (*node)(0);
        const Number exact = left_value + (right_value - left_value) * x;
        LIBMESH_ASSERT_NUMBERS_EQUAL(
            exact, (*sys.solution)(node->dof_number(sys.number(), u, 0)), TOLERANCE * TOLERANCE);
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(HeterogeneousNonlinearConstraintTest);
