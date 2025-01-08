#include <libmesh/libmesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/system.h>

#include <timpi/parallel_implementation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class ConnectedComponentsTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to ensure that the count returned by
   * MeshTools::n_connected_components() is correct for hand-checked
   * cases.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ConnectedComponentsTest );

  /**
   * These tests require a valid solver package to be enabled because
   * they build a SparseMatrix.
   * Note: LIBMESH_HAVE_SOLVER is actually defined in libmesh/system.h
   */
#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testEdge2 );
  CPPUNIT_TEST( testEdge3 );
  CPPUNIT_TEST( testEdge4 );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  void testEdge(unsigned n_elem,
                ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    // We don't yet support n_connected_components() on distributed
    // meshes
    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_line(mesh,
                                      n_elem,
                                      /*xmin=*/0.,
                                      /*xmax=*/1.,
                                      elem_type);

    const dof_id_type n_components_orig =
      MeshTools::n_connected_components(mesh);

    CPPUNIT_ASSERT_EQUAL(n_components_orig, dof_id_type(1));

    // We can't disconnect a too-small mesh
    if (n_elem < 3)
      return;

    for (const auto & elem : mesh.element_ptr_range())
      {
        const Point c = elem->vertex_average();

        if (std::abs(c(0) - 0.5) < 0.2)
          mesh.delete_elem(elem);
      }

    mesh.prepare_for_use();

    const dof_id_type n_components_broken =
      MeshTools::n_connected_components(mesh);

    CPPUNIT_ASSERT_EQUAL(n_components_broken, dof_id_type(2));

    // Add a constraint connecting the original end nodes
    auto matrix = SparseMatrix<Number>::build (mesh.comm());

    dof_id_type left_node = 0, right_node = 0;
    for (const auto & elem : mesh.element_ptr_range())
      for (const auto & node : elem->node_ref_range())
        {
          if (node(0) == Real(0))
            left_node = node.id();
          if (node(0) == Real(1))
            right_node = node.id();
        }

    mesh.comm().max(left_node);
    mesh.comm().max(right_node);

    const dof_id_type n_nodes = mesh.n_nodes();
    CPPUNIT_ASSERT_EQUAL(n_nodes, mesh.max_node_id());

    processor_id_type my_rank = mesh.comm().rank();
    matrix->init(n_nodes, n_nodes-1,
                 (my_rank ? 0 : n_nodes),
                 (my_rank ? 0 : n_nodes-1),
                 (my_rank ? 0 : n_nodes-1),
                 0);

    if (!my_rank)
      {
        for (auto i : make_range(n_nodes))
          {
            if (i < right_node)
              matrix->set(i,i,1);
            if (i > right_node)
              matrix->set(i,i-1,1);
          }
        if (left_node > right_node)
          matrix->set(right_node, left_node-1, 1);
        else
          matrix->set(right_node, left_node, 1);
      }

    matrix->close();

    mesh.copy_constraint_rows(*matrix);

    const dof_id_type n_components_constrained =
      MeshTools::n_connected_components(mesh);

    CPPUNIT_ASSERT_EQUAL(n_components_constrained, dof_id_type(1));
  }

public:
  void setUp() {}

  void tearDown() {}

  void testEdge2()
  {
    LOG_UNIT_TEST;

    testEdge(/*n_elem=*/10, EDGE2);
  }

  void testEdge3()
  {
    LOG_UNIT_TEST;

    testEdge(/*n_elem=*/10, EDGE3);
  }

  void testEdge4()
  {
    LOG_UNIT_TEST;

    testEdge(/*n_elem=*/10, EDGE4);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( ConnectedComponentsTest );
