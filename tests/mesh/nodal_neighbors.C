#include <libmesh/libmesh.h>
#include <libmesh/node.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/utility.h>
#include <libmesh/reference_elem.h>
#include <unordered_map>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class NodalNeighborsTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to ensure that MeshTools::find_nodal_neighbors()
   * works as expected.  If the numbering of MeshGeneration::build_line()
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

  CPPUNIT_TEST( testTri3 );
  CPPUNIT_TEST( testTri6 );
  CPPUNIT_TEST( testTri7 );

  CPPUNIT_TEST( testQuad4 );
  CPPUNIT_TEST( testQuad8 );
  CPPUNIT_TEST( testQuad9 );

  CPPUNIT_TEST( testTet4 );
  CPPUNIT_TEST( testTet10 );
  CPPUNIT_TEST( testTet14 );

  CPPUNIT_TEST( testPyramid5 );
  CPPUNIT_TEST( testPyramid13 );
  CPPUNIT_TEST( testPyramid14 );
  CPPUNIT_TEST( testPyramid18 );

  CPPUNIT_TEST( testPrism6 );
  CPPUNIT_TEST( testPrism18 );
  CPPUNIT_TEST( testPrism20 );
  CPPUNIT_TEST( testPrism21 );

  CPPUNIT_TEST( testHex8 );
  CPPUNIT_TEST( testHex20 );
  CPPUNIT_TEST( testHex27 );

  CPPUNIT_TEST_SUITE_END();

private:

  // Hard-coded neighbor relationships for each element type
  using Map = std::map<ElemType, std::map<dof_id_type, std::set<dof_id_type>>>;

  static Map build_elem_type_to_neighbor_map()
  {
    Map m;
    //m.reserve(10);

    m[EDGE2][0].insert(1);
    m[EDGE2][1].insert(0);

    m[EDGE3][0].insert(2);
    m[EDGE3][1].insert(2);
    m[EDGE3][2].insert({0, 1});

    m[EDGE4][0].insert(2);
    m[EDGE4][1].insert(3);
    m[EDGE4][2].insert({0, 3});
    m[EDGE4][3].insert({2, 1});

    m[TRI3][0].insert({1, 2});
    m[TRI3][1].insert({0, 2});
    m[TRI3][2].insert({1, 0});

    m[TRI6][0].insert({3, 5});
    m[TRI6][1].insert({3, 4});
    m[TRI6][2].insert({5, 4});
    m[TRI6][3].insert({0, 1});
    m[TRI6][4].insert({1, 2});
    m[TRI6][5].insert({2, 0});

    m[TRI7] = m[TRI6];
    // Node 6 of TRI7 is a face node, has no neighbors because it is not on an edge

    m[QUAD4][0].insert({1, 3});
    m[QUAD4][1].insert({0, 2});
    m[QUAD4][2].insert({1, 3});
    m[QUAD4][3].insert({0, 2});

    m[QUAD8][0].insert({4, 7});
    m[QUAD8][1].insert({4, 5});
    m[QUAD8][2].insert({5, 6});
    m[QUAD8][3].insert({6, 7});
    m[QUAD8][4].insert({0, 1});
    m[QUAD8][5].insert({1, 2});
    m[QUAD8][6].insert({2, 3});
    m[QUAD8][7].insert({0, 3});

    m[QUAD9] = m[QUAD8];
    // Node 8 of QUAD9 is a face node, has no neighbors because it is not on an edge

    m[TET4][0].insert({1, 2, 3});
    m[TET4][1].insert({0, 2, 3});
    m[TET4][2].insert({0, 1, 3});
    m[TET4][3].insert({0, 1, 2});

    m[TET10][0].insert({4, 6, 7});
    m[TET10][1].insert({4, 5, 8});
    m[TET10][2].insert({5, 6, 9});
    m[TET10][3].insert({7, 8, 9});
    m[TET10][4].insert({0, 1});
    m[TET10][5].insert({1, 2});
    m[TET10][6].insert({0, 2});
    m[TET10][7].insert({0, 3});
    m[TET10][8].insert({1, 3});
    m[TET10][9].insert({2, 3});

    m[TET14] = m[TET10];
    // Nodes 10-13 are face nodes and have no neighbors

    m[PYRAMID5][0].insert({1, 3, 4});
    m[PYRAMID5][1].insert({0, 2, 4});
    m[PYRAMID5][2].insert({1, 3, 4});
    m[PYRAMID5][3].insert({0, 2, 4});
    m[PYRAMID5][4].insert({0, 1, 2, 3});

    m[PYRAMID13][0].insert({5, 8, 9});
    m[PYRAMID13][1].insert({5, 6, 10});
    m[PYRAMID13][2].insert({6, 7, 11});
    m[PYRAMID13][3].insert({7, 8, 12});
    m[PYRAMID13][4].insert({9, 10, 11, 12});
    m[PYRAMID13][5].insert({0, 1});
    m[PYRAMID13][6].insert({1, 2});
    m[PYRAMID13][7].insert({2, 3});
    m[PYRAMID13][8].insert({0, 3});
    m[PYRAMID13][9].insert({0, 4});
    m[PYRAMID13][10].insert({1, 4});
    m[PYRAMID13][11].insert({2, 4});
    m[PYRAMID13][12].insert({3, 4});

    m[PYRAMID14] = m[PYRAMID13];
    // Node 13 is a face node and has no neighbors

    m[PYRAMID18] = m[PYRAMID14];
    // Nodes 14-17 are face nodes and have no neighbors

    m[PRISM6][0].insert({1, 2, 3});
    m[PRISM6][1].insert({0, 2, 4});
    m[PRISM6][2].insert({0, 1, 5});
    m[PRISM6][3].insert({0, 4, 5});
    m[PRISM6][4].insert({1, 3, 5});
    m[PRISM6][5].insert({2, 3, 4});

    m[PRISM18][0].insert({6, 8, 9});
    m[PRISM18][1].insert({6, 7, 10});
    m[PRISM18][2].insert({7, 8, 11});
    m[PRISM18][3].insert({9, 12, 14});
    m[PRISM18][4].insert({10, 12, 13});
    m[PRISM18][5].insert({11, 13, 14});
    m[PRISM18][6].insert({0, 1});
    m[PRISM18][7].insert({1, 2});
    m[PRISM18][8].insert({0, 2});
    m[PRISM18][9].insert({0, 3});
    m[PRISM18][10].insert({1, 4});
    m[PRISM18][11].insert({2, 5});
    m[PRISM18][12].insert({3, 4});
    m[PRISM18][13].insert({4, 5});
    m[PRISM18][14].insert({3, 5});
    // Nodes 15-17 are face nodes and have no neighbors

    m[PRISM20] = m[PRISM18];
    // Nodes 18-19 are face nodes and have no neighbors

    m[PRISM21] = m[PRISM20];
    // Node 20 is an interior node and has no neighbors

    m[HEX8][0].insert({1, 3, 4});
    m[HEX8][1].insert({0, 2, 5});
    m[HEX8][2].insert({1, 3, 6});
    m[HEX8][3].insert({0, 2, 7});
    m[HEX8][4].insert({0, 5, 7});
    m[HEX8][5].insert({1, 4, 6});
    m[HEX8][6].insert({2, 5, 7});
    m[HEX8][7].insert({3, 4, 6});

    m[HEX20][0].insert({8, 11, 12});
    m[HEX20][1].insert({8, 9, 13});
    m[HEX20][2].insert({9, 10, 14});
    m[HEX20][3].insert({10, 11, 15});
    m[HEX20][4].insert({12, 16, 19});
    m[HEX20][5].insert({13, 16, 17});
    m[HEX20][6].insert({14, 17, 18});
    m[HEX20][7].insert({15, 18, 19});
    m[HEX20][8].insert({0, 1});
    m[HEX20][9].insert({1, 2});
    m[HEX20][10].insert({2, 3});
    m[HEX20][11].insert({0, 3});
    m[HEX20][12].insert({0, 4});
    m[HEX20][13].insert({1, 5});
    m[HEX20][14].insert({2, 6});
    m[HEX20][15].insert({3, 7});
    m[HEX20][16].insert({4, 5});
    m[HEX20][17].insert({5, 6});
    m[HEX20][18].insert({6, 7});
    m[HEX20][19].insert({4, 7});

    m[HEX27] = m[HEX20];
    // Nodes 20-26 are face nodes and have no neighbors

    return m;
  };

  inline static const Map elem_type_to_neighbor_map = build_elem_type_to_neighbor_map();

protected:

  // Builds a 1D mesh with the specified ElemType and number of elements
  void do_test(ElemType elem_type)
  {
    const auto * ref_elem = &(ReferenceElem::get(elem_type));
    const auto dim = ref_elem->dim();

    ReplicatedMesh mesh(*TestCommWorld, dim);

    unsigned int n_elems_per_side = 4 / Elem::type_to_default_order_map[elem_type];

    switch (dim)
      {
        case 1:
          MeshTools::Generation::build_line(mesh, n_elems_per_side, 0., 1., elem_type);
          break;
        case 2:
          MeshTools::Generation::build_square(
              mesh, n_elems_per_side, n_elems_per_side, 0., 1., 0., 1., elem_type);
          break;

        case 3:
          MeshTools::Generation::build_cube(mesh,
                                            n_elems_per_side,
                                            n_elems_per_side,
                                            n_elems_per_side,
                                            0.,
                                            1.,
                                            0.,
                                            1.,
                                            0.,
                                            1.,
                                            elem_type);
          break;

        default:
          libmesh_error_msg("Unsupported dimension " << dim);
      }

    const auto & neighbor_map = libmesh_map_find(elem_type_to_neighbor_map, elem_type);

    // find_nodal_neighbors() needs a data structure which is prepared by another function
    std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
    MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

    // Loop over the nodes and call find_nodal_neighbors()
    std::vector<const Node*> neighbor_nodes;

    std::set<dof_id_type> node_ids_checked;
    for (const auto * elem : mesh.element_ptr_range())
      {
        for (const auto & node : elem->node_ref_range())
          {
            // Have we already checked this node?
            if (node_ids_checked.find(node.id()) != node_ids_checked.end())
              continue;

            MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbor_nodes);

            for (const auto * neigh : neighbor_nodes)
              {
                // Find the common elements between node and neigh
                const auto & elems_containing_node = libmesh_map_find(nodes_to_elem_map, node.id());
                const auto & elems_containing_neigh = libmesh_map_find(nodes_to_elem_map, neigh->id());
                const Elem * common_elem = nullptr;
                for (const auto * neigh_elem : elems_containing_neigh)
                  {
                    if ((std::find(elems_containing_node.begin(), elems_containing_node.end(), neigh_elem) !=
                         elems_containing_node.end()))
                      common_elem = neigh_elem;
                    else
                      continue;

                    // Which local elem node numbers are node and neighbor?
                    unsigned node_id_local = common_elem->local_node(node.id());
                    unsigned neigh_id_local = common_elem->local_node(neigh->id());

                    // For a given element type, we know the hard-coded neighbor relationship
                    // between local node ids. Make sure the neigbor from find_nodal_neighbors
                    // was intended to be a neighbor. Note that this test does not currently check
                    // for neighbors that find_nodal_neighbors may not have found.
                    const auto & correct_neighbor_ids = libmesh_map_find(neighbor_map, node_id_local);
                    CPPUNIT_ASSERT(correct_neighbor_ids.find(neigh_id_local) != correct_neighbor_ids.end());

                  }
              }

            node_ids_checked.insert(node.id());
          }
      }
  }

public:
  void setUp() {}

  void tearDown() {}

  void testEdge2()
  {
    LOG_UNIT_TEST;
    do_test(EDGE2);
  }


  void testEdge3()
  {
    LOG_UNIT_TEST;
    do_test(EDGE3);
  }


  void testEdge4()
  {
    LOG_UNIT_TEST;
    do_test(EDGE4);
  }

  void testTri3()
  {
    LOG_UNIT_TEST;
    do_test(TRI3);
  }

  void testTri6()
  {
    LOG_UNIT_TEST;
    do_test(TRI6);
  }

  void testTri7()
  {
    LOG_UNIT_TEST;
    do_test(TRI7);
  }

  void testQuad4()
  {
    LOG_UNIT_TEST;
    do_test(QUAD4);
  }

  void testQuad8()
  {
    LOG_UNIT_TEST;
    do_test(QUAD8);
  }

  void testQuad9()
  {
    LOG_UNIT_TEST;
    do_test(QUAD9);
  }

  void testTet4()
  {
    LOG_UNIT_TEST;
    do_test(TET4);
  }

  void testTet10()
  {
    LOG_UNIT_TEST;
    do_test(TET10);
  }

  void testTet14()
  {
    LOG_UNIT_TEST;
    do_test(TET14);
  }

  void testPyramid5()
  {
    LOG_UNIT_TEST;
    do_test(PYRAMID5);
  }

  void testPyramid13()
  {
    LOG_UNIT_TEST;
    do_test(PYRAMID13);
  }

  void testPyramid14()
  {
    LOG_UNIT_TEST;
    do_test(PYRAMID14);
  }

  void testPyramid18()
  {
    LOG_UNIT_TEST;
    do_test(PYRAMID18);
  }

  void testPrism6()
  {
    LOG_UNIT_TEST;
    do_test(PRISM6);
  }

  void testPrism18()
  {
    LOG_UNIT_TEST;
    do_test(PRISM18);
  }

  void testPrism20()
  {
    LOG_UNIT_TEST;
    do_test(PRISM20);
  }

  void testPrism21()
  {
    LOG_UNIT_TEST;
    do_test(PRISM21);
  }

  void testHex8()
  {
    LOG_UNIT_TEST;
    do_test(HEX8);
  }

  void testHex20()
  {
    LOG_UNIT_TEST;
    do_test(HEX20);
  }

  void testHex27()
  {
    LOG_UNIT_TEST;
    do_test(HEX27);
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
        elem->set_node(0, mesh.node_ptr(conn[2*e] - 1));     // convert to 0-based index
        elem->set_node(1, mesh.node_ptr(conn[2*e + 1] - 1)); // convert to 0-based index
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
