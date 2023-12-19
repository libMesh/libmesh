#include <libmesh/libmesh.h>

#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include <timpi/parallel_implementation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MixedOrderTest : public CppUnit::TestCase
{
  /**
   * The goal of this test set is to test libMesh manipulation of
   * meshes with mixed Elem orders.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MixedOrderTest );

  CPPUNIT_TEST( testFindNeighbors );
  CPPUNIT_TEST( testStitch );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  static int n_neighbor_links (const MeshBase & mesh)
  {
    int n_neighbors = 0;
    for (const auto & elem : mesh.local_element_ptr_range())
      for (const auto & neigh : elem->neighbor_ptr_range())
        if (neigh)
          ++n_neighbors;
    mesh.comm().sum(n_neighbors);
    return n_neighbors;
  }


  void testFindNeighbors()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    // Construct a multi-element Tri3 mesh
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/3, /*ny=*/3,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/TRI3);

    // Make some of them into Tri6
    auto range_start = mesh.elements_begin();
    const auto range_end = mesh.elements_end();
    while (range_start != range_end && (*range_start)->id() < 5)
      ++range_start;

    const int old_n_neighbors = n_neighbor_links(mesh);

    mesh.all_second_order_range({range_start,range_end},
                                /*full_ordered=*/true);
    const int new_n_neighbors = n_neighbor_links(mesh);
    CPPUNIT_ASSERT_EQUAL(old_n_neighbors, new_n_neighbors);

    mesh.find_neighbors();
    const int newer_n_neighbors = n_neighbor_links(mesh);
    CPPUNIT_ASSERT_EQUAL(old_n_neighbors, newer_n_neighbors);
  }

  void testStitch()
  {
    LOG_UNIT_TEST;

    Mesh mesh0(*TestCommWorld),
         mesh1(*TestCommWorld);

    // Quad4 mesh and Quad9 mesh
    MeshTools::Generation::build_square(mesh0,
                                        /*nx=*/3, /*ny=*/3,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/QUAD4);
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), dof_id_type(16));

    MeshTools::Generation::build_square(mesh1,
                                        /*nx=*/3, /*ny=*/3,
                                        /*xmin=*/-1., /*xmax=*/0.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/QUAD9);
    CPPUNIT_ASSERT_EQUAL(mesh1.n_nodes(), dof_id_type(49));

    // Stitch them together
    mesh0.stitch_meshes(mesh1, 3, 1, TOLERANCE, /*clear_bcids*/ true,
                        /*verbose*/ false, /*binary_search*/ false,
                        /*enforce_all_nodes_match*/ false,
                        /*merge_all_or_nothing*/ false);
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), dof_id_type(61));

    const int stitched_n_neighbors = n_neighbor_links(mesh0);
    CPPUNIT_ASSERT_EQUAL(stitched_n_neighbors, 54);

    mesh0.all_second_order_range(mesh0.element_ptr_range(),
                                 /*full_ordered=*/true);
    const int new_n_neighbors = n_neighbor_links(mesh0);
    CPPUNIT_ASSERT_EQUAL(new_n_neighbors, 54);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MixedOrderTest );
