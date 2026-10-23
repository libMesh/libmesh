#include <libmesh/constraint_elimination_dof_map.h>
#include <libmesh/const_function.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/reduced_sparse_matrix.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
/**
 * Assembles a non-symmetric operator, a Laplacian plus an advection term, so that the test can
 * tell a constrained row of C^T K C from a constrained column
 */
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

#if defined(LIBMESH_HAVE_PETSC) && defined(LIBMESH_ENABLE_AMR) &&                                  \
    defined(LIBMESH_ENABLE_PERIODIC) && defined(LIBMESH_ENABLE_DIRICHLET) && LIBMESH_DIM > 1
  CPPUNIT_TEST(testReducedOperator);
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
};

CPPUNIT_TEST_SUITE_REGISTRATION(ConstraintEliminationTest);
