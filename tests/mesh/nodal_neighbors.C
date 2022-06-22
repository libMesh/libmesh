#include <libmesh/libmesh.h>
#include <libmesh/node.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class NodalNeighborsTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to ensure that MeshTools::find_nodal_neighbors()
   * works in 1D.  If the numbering of MeshGeneration::build_line()
   * ever changes, this test will break, as it compares hand-checked
   * hard-coded "validation" data with the results of
   * MeshTools::find_nodal_neighbors().  We also use a ReplicatedMesh
   * here to match the hard-coded numbering.
   *
   * The testOrientation() test is not specifically related to
   * find_nodal_neighbors(), instead it is checking that we can still
   * find_neighbors() correctly in 1D when the mesh is topologically a
   * straight line, but not all elements have the same "orientation"
   * (as defined by their local node numbering). As far as I know, we
   * don't require 2D/3D elements to have the same orientation in
   * order for them to be considered neighbors, so this test ensures
   * the same thing works for 1D elements.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( NodalNeighborsTest );

  CPPUNIT_TEST( testEdge2 );
  CPPUNIT_TEST( testEdge3 );
  CPPUNIT_TEST( testEdge4 );
  CPPUNIT_TEST( testOrientation );

  CPPUNIT_TEST_SUITE_END();

protected:

  // Builds a 1D mesh with the specified ElemType and number of elements
  void do_test(unsigned n_elem,
               ElemType elem_type,
               dof_id_type * validation_data)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/1);

    MeshTools::Generation::build_line(mesh,
                                      n_elem,
                                      /*xmin=*/0.,
                                      /*xmax=*/1.,
                                      elem_type);

    // find_nodal_neighbors() needs a data structure which is prepared by another function
    std::vector<std::vector<const Elem *>> nodes_to_elem_map;
    MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

    // Loop over the nodes and call find_nodal_neighbors()
    {
      std::vector<const Node*> neighbor_nodes;

      unsigned ctr = 0;
      for (const auto & node : mesh.node_ptr_range())
        {
          MeshTools::find_nodal_neighbors(mesh, *node, nodes_to_elem_map, neighbor_nodes);

          // The entries in neighbor_nodes are just sorted according
          // to memory address, which is somewhat arbitrary, so create
          // a vector sorted by IDs for test purposes.
          std::vector<dof_id_type> neighbor_node_ids(neighbor_nodes.size());
          for (std::size_t i=0; i<neighbor_nodes.size(); ++i)
            neighbor_node_ids[i] = neighbor_nodes[i]->id();
          std::sort(neighbor_node_ids.begin(), neighbor_node_ids.end());

          // Compare to validation_data
          for (std::size_t j=0; j<neighbor_node_ids.size(); ++j)
            {
              CPPUNIT_ASSERT_EQUAL( validation_data[2*ctr + j], neighbor_node_ids[j] );
            }

          ++ctr;
        }
    }
  }

public:
  void setUp() {}

  void tearDown() {}

  void testEdge2()
  {
    LOG_UNIT_TEST;

    // 11 nodes, 2 neighbor entries per node
    dof_id_type validation_data[22] =
      {
        1, DofObject::invalid_id,
        0, 2,
        1, 3,
        2, 4,
        3, 5,
        4, 6,
        5, 7,
        6, 8,
        7, 9,
        8, 10,
        9, DofObject::invalid_id
      };

    do_test(/*n_elem=*/10, EDGE2, validation_data);
  }


  void testEdge3()
  {
    LOG_UNIT_TEST;

    // 11 nodes, 2 neighbor entries per node
    dof_id_type validation_data[22] =
      {
        2, DofObject::invalid_id,
        2, 4,
        0, 1,
        4, 6,
        1, 3,
        6, 8,
        3, 5,
        8, 10,
        5, 7,
        10, DofObject::invalid_id,
        7, 9
      };

    do_test(/*n_elem=*/5, EDGE3, validation_data);
  }


  void testEdge4()
  {
    LOG_UNIT_TEST;

    // 10 nodes, 2 neighbor entries per node
    dof_id_type validation_data[20] =
      {
        2, DofObject::invalid_id,
        3, 5,
        0, 3,
        1, 2,
        6, 8,
        1, 6,
        4, 5,
        9, DofObject::invalid_id,
        4, 9,
        7, 8
      };

    do_test(/*n_elem=*/3, EDGE4, validation_data);
  }

  void testOrientation()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);

    // The (1-based) node numbering and element orientation (represented
    // by arrows) for this mesh:
    // 1 -> 3 -> 4 -> 7 -> 8 <- 6 <- 5 <- 2
    // So the two elements that meet at node 8 have opposite orientation
    // and were not detected as neighbors using the original (before the
    // addition of this test) find_neighbors() algorithm.

    // These nodes are copied from an exo file where we originally
    // noticed the problem, but otherwise are not significant to the
    // test. Note: the actual ids are 0-based but we are using a
    // 1-based connectivity array below which is copied from the exo
    // file.
    mesh.add_point(Point(1.68, -1.695, 8.298), /*id=*/0);
    mesh.add_point(Point(5.55, -1.695, 8.298), /*id=*/1);
    mesh.add_point(Point(1.68, -0.175, 8.298), /*id=*/2);
    mesh.add_point(Point(1.68, -0.175, 9.643), /*id=*/3);
    mesh.add_point(Point(5.55, -0.175, 8.298), /*id=*/4);
    mesh.add_point(Point(5.55, -0.175, 9.643), /*id=*/5);
    mesh.add_point(Point(1.68, -0.075, 9.643), /*id=*/6);
    mesh.add_point(Point(5.55, -0.075, 9.643), /*id=*/7);

    // 1-based connectivity array (2 nodes per Elem) copied directly
    // from exo file.  We will convert these to 0-based ids when they
    // are used.
    std::vector<unsigned int> conn =
      {
        7, 8,
        1, 3,
        2, 5,
        3, 4,
        5, 6,
        4, 7,
        6, 8
      };

    // Add 7 EDGE2 elements and assign connectivity
    for (unsigned int e=0; e<7; ++e)
      {
        Elem * elem = mesh.add_elem(Elem::build_with_id(EDGE2, e));
        elem->set_node(0) = mesh.node_ptr(conn[2*e] - 1);     // convert to 0-based index
        elem->set_node(1) = mesh.node_ptr(conn[2*e + 1] - 1); // convert to 0-based index
      }

    // Find neighbors, etc.
    mesh.prepare_for_use();

    for (const auto & elem : mesh.element_ptr_range())
      {
        auto elem_id = elem->id();

        // Elems 1, 2 should have no neighbor on side 0
        if (elem_id == 1 || elem_id == 2)
          {
            CPPUNIT_ASSERT(elem->neighbor_ptr(0) == nullptr);
            CPPUNIT_ASSERT(elem->neighbor_ptr(1) != nullptr);
          }
        // Otherwise, elem should have neighbor on both sides
        else
          {
            CPPUNIT_ASSERT(elem->neighbor_ptr(0) != nullptr);
            CPPUNIT_ASSERT(elem->neighbor_ptr(1) != nullptr);
          }

        // Debugging
        // libMesh::out << "elem " << elem_id << std::endl;
        // for (auto side_idx : make_range(elem->n_sides()))
        //   {
        //     if (elem->neighbor_ptr(side_idx))
        //       libMesh::out << "Neighbor on side " << side_idx << std::endl;
        //     else
        //       libMesh::out << "No neighbor on side " << side_idx << std::endl;
        //   }
      }
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( NodalNeighborsTest );
