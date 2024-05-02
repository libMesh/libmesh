#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_norm_type.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exact_solution.h>
#include <libmesh/int_range.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/parallel.h>
#include <libmesh/quadrature.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/petsc_macro.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <string>

using namespace libMesh;

namespace {

// A simple u"+u=f
void assemble_matrix_and_rhs(EquationSystems& es,
                             const std::string&)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");
  const DofMap& dof_map = system.get_dof_map();
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  FEMContext c(system);

  const unsigned short dim = mesh.mesh_dimension();
  FEBase * fe = nullptr;
  c.get_element_fe(0, fe, dim);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<Point> & xyz = fe->get_xyz();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  std::vector<dof_id_type> & dof_indices = c.get_dof_indices();

  // The subvectors and submatrices we need to fill:
  DenseMatrix<Number> & K = c.get_elem_jacobian();
  DenseVector<Number> & F = c.get_elem_residual();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      if(elem->type() == NODEELEM)
        continue;

      c.pre_fe_reinit(system, elem);
      c.elem_fe_reinit();

      unsigned int n_dofs = dof_indices.size();
      unsigned int n_qpoints = c.get_element_qrule().n_points();

      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          const Gradient grad_u = c.interior_gradient(0, qp);
          const Number u = c.interior_value(0, qp);
          const Point p = xyz[qp];

          const Number forcing = p(0);

          for (unsigned int i=0; i != n_dofs; i++)
            {
              F(i) += JxW[qp] * (grad_u * dphi[i][qp] + (forcing-u) * phi[i][qp]);
              for (unsigned int j=0; j != n_dofs; ++j)
                K(i,j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp] - phi[i][qp]*phi[j][qp]);
            }
        }

      dof_map.constrain_element_matrix_and_vector (K, F, dof_indices);

      matrix.add_matrix (K, dof_indices);
      system.rhs->add_vector (F, dof_indices);
    }

  system.rhs->close();
  matrix.close();
}
}



class ConstraintOperatorTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ConstraintOperatorTest );

#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( test1DCoarseningOperator );
  CPPUNIT_TEST( test1DCoarseningNewNodes );

#ifdef LIBMESH_HAVE_EXODUS_API
  // These tests seem to fail on older PETScs with the following error:
  // [0]PETSC ERROR: Nonconforming object sizes
  // [0]PETSC ERROR: Square matrix only for in-place
  // [0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
  // [0]PETSC ERROR: Petsc Release Version 3.8.0, unknown
  // So we only run them if you have new enough PETSc. We
  // don't know exactly how new of a PETSc is required, though,
  // so for now it's anything greater than 3.8. See also:
  // libMesh/libmesh#3803.
#if PETSC_RELEASE_GREATER_THAN(3,8,0)
  CPPUNIT_TEST( testCoreformGeom2 );
  CPPUNIT_TEST( testCoreformGeom4 );
  CPPUNIT_TEST( testCoreformGeom8 );
#endif

#endif
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}


  void test1DCoarseningOperator ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);

    ExplicitSystem &sys =
      es.add_system<LinearImplicitSystem> ("test");

    sys.add_variable("u", FIRST);
    sys.attach_assemble_function (assemble_matrix_and_rhs);

    // Make sure numbering matches our constraint operator matrix
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_line (mesh,
                                       10,
                                       0., 1.,
                                       EDGE2);

    // We should be prepared at this point
    CPPUNIT_ASSERT(mesh.is_prepared());
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    auto matrix = SparseMatrix<Number>::build (mesh.comm());

    // Create a projection matrix from 10 elements to 5.  We might as
    // well just make it on processor 0.
    processor_id_type my_rank = mesh.comm().rank();
    matrix->init(11, 6,
                 (my_rank ? 0 : 11),
                 (my_rank ? 0 : 6),
                 (my_rank ? 0 : 6),
                 0);

    if (!my_rank)
      {
        for (auto i : make_range(5))
          {
            matrix->set(i*2, i, 1);
            matrix->set(i*2+1, i, 0.5);
            matrix->set(i*2+1, i+1, 0.5);
          }
        matrix->set(10, 5, 1);
      }
    matrix->close();

    mesh.copy_constraint_rows(*matrix);

    // We should still be prepared at this point
    CPPUNIT_ASSERT(mesh.is_prepared());
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    // Do a redundant prepare-for-use to catch any issues there
    mesh.prepare_for_use();

    es.init ();
    sys.solve();

    const DofMap & dof_map = sys.get_dof_map();

    CPPUNIT_ASSERT_EQUAL(dof_id_type(5), dof_map.n_constrained_dofs());

    // Now compare the results to just solving on 5 elements.

    // Matching up partitioning is hard, so just replicate+serialize
    ReplicatedMesh mesh2(*TestCommWorld);
    EquationSystems es2(mesh);

    ExplicitSystem &sys2 =
      es2.add_system<LinearImplicitSystem> ("test");

    sys2.add_variable("u", FIRST);
    sys2.attach_assemble_function (assemble_matrix_and_rhs);

    MeshTools::Generation::build_line (mesh2,
                                       5,
                                       0., 1.,
                                       EDGE2);

    es2.init ();
    sys2.solve();

    // Integrate any differences on the fine grid
    ExactSolution exact(es);

    // Serialize rather than figure out proper ghosting for an
    // arbitrary partitioning mix
    auto serialized_solution2 =
      NumericVector<Number>::build(*TestCommWorld);
    serialized_solution2->init(sys2.solution->size(), false, SERIAL);
    sys2.solution->localize(*serialized_solution2);

    MeshFunction coarse_solution
      (es2, *serialized_solution2, sys2.get_dof_map(), 0);

    exact.attach_exact_value(0, &coarse_solution);
    exact.compute_error("test", "u");
    Real err = exact.l2_error("test", "u");
    CPPUNIT_ASSERT_LESS(Real(TOLERANCE*TOLERANCE), err);
  }


  void testCoreform(const std::string & meshfile,
                    const std::string & matrixfile,
                    Real expected_norm,
                    Order order=FIRST)
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);

    ExplicitSystem &sys =
      es.add_system<LinearImplicitSystem> ("test");

    sys.add_variable("u", order);
    sys.attach_assemble_function (assemble_matrix_and_rhs);

    // Make sure numbering matches our constraint operator matrix
    mesh.allow_renumbering(false);

    mesh.read(meshfile);

    CPPUNIT_ASSERT(mesh.is_prepared());

    // For these matrices Coreform has been experimenting with PETSc
    // solvers which take the transpose of what we expect, so we'll
    // un-transpose here.
    auto matrix = SparseMatrix<Number>::build (mesh.comm());
    matrix->read_matlab(matrixfile);
    matrix->get_transpose(*matrix);

    mesh.copy_constraint_rows(*matrix);

    // We should be prepared at this point
    CPPUNIT_ASSERT(mesh.is_prepared());
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    // Do a redundant prepare-for-use to catch any issues there
    mesh.prepare_for_use();

    es.init ();
    sys.solve();

    // Skip norm calculations on NodeElems
    std::set<unsigned int> skip_dimensions {0};
    Real h1_norm = sys.calculate_norm(*sys.solution, 0, H1,
                                      &skip_dimensions);

    LIBMESH_ASSERT_FP_EQUAL(h1_norm, expected_norm, TOLERANCE*std::sqrt(TOLERANCE));
  }


  void testCoreformGeom2()
  {
    LOG_UNIT_TEST;
    testCoreform("meshes/geom_2.exo", "matrices/geom_2_extraction_op.m", 0.0959827039398013);
  }


  void testCoreformGeom4()
  {
    LOG_UNIT_TEST;
    testCoreform("meshes/geom_4.exo", "matrices/geom_4_extraction_op.m", 0.10404051279393);
  }


  void testCoreformGeom8()
  {
    LOG_UNIT_TEST;
    testCoreform("meshes/geom_8.exo", "matrices/geom_8_extraction_op.m", 0.103252835327016);
  }


  void test1DCoarseningNewNodes ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);

    ExplicitSystem &sys =
      es.add_system<LinearImplicitSystem> ("test");

    sys.add_variable("u", FIRST);
    sys.attach_assemble_function (assemble_matrix_and_rhs);

    // Make sure numbering matches our constraint operator matrix
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_line (mesh,
                                       10,
                                       0., 1.,
                                       EDGE2);

    // We should be prepared at this point
    CPPUNIT_ASSERT(mesh.is_prepared());
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    auto matrix = SparseMatrix<Number>::build (mesh.comm());

    // Create a projection matrix from 10 elements to 4.  We'll only
    // be leaving the end nodes and one central node unconstrained, so
    // we should see two new unconstrained NodeElem added.
    matrix->read_matlab("meshes/constrain10to4.m");

    mesh.copy_constraint_rows(*matrix);

    // We should still be prepared at this point
    CPPUNIT_ASSERT(mesh.is_prepared());
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    // Do a redundant prepare-for-use to catch any issues there
    mesh.prepare_for_use();

    es.init ();
    sys.solve();

    const DofMap & dof_map = sys.get_dof_map();

    CPPUNIT_ASSERT_EQUAL(dof_id_type(11+2), mesh.n_nodes());
    CPPUNIT_ASSERT_EQUAL(dof_id_type(10+2), mesh.n_elem());
    CPPUNIT_ASSERT_EQUAL(dof_id_type(8), dof_map.n_constrained_dofs());
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ConstraintOperatorTest );
