// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh.h>
#include <libmesh/serial_mesh.h>
#include <libmesh/node.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>

#include "test_comm.h"

using namespace libMesh;

class NodalNeighborsTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to ensure that MeshTools::find_nodal_neighbors()
   * works in 1D.  If the numbering of MeshGeneration::build_line()
   * ever changes, this test will break, as it compares hand-checked
   * hard-coded "validation" data with the results of
   * MeshTools::find_nodal_neighbors().  We also use a SerialMesh here
   * to match the hard-coded numbering.
   */
public:
  CPPUNIT_TEST_SUITE( NodalNeighborsTest );

  CPPUNIT_TEST( testEdge2 );
  CPPUNIT_TEST( testEdge3 );
  CPPUNIT_TEST( testEdge4 );

  CPPUNIT_TEST_SUITE_END();

protected:

  // Builds a 1D mesh with the specified ElemType and number of elements
  void do_test(unsigned n_elem,
               ElemType elem_type,
               dof_id_type * validation_data)
  {
    SerialMesh mesh(*TestCommWorld, /*dim=*/1);

    MeshTools::Generation::build_line(mesh,
                                      n_elem,
                                      /*xmin=*/0.,
                                      /*xmax=*/1.,
                                      elem_type);

    // find_nodal_neighbors() needs a data structure which is prepared by another function
    std::vector<std::vector<const Elem*> > nodes_to_elem_map;
    MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

    // Loop over the nodes and call find_nodal_neighbors()
    {
      MeshBase::const_node_iterator       nd     = mesh.nodes_begin();
      const MeshBase::const_node_iterator end_nd = mesh.nodes_end();

      std::vector<const Node*> neighbor_nodes;

      unsigned ctr = 0;
      for (; nd != end_nd; ++nd, ++ctr)
        {
          Node* node = *nd;

          MeshTools::find_nodal_neighbors(mesh, *node, nodes_to_elem_map, neighbor_nodes);

          // The entries in neighbor_nodes are just sorted according
          // to memory address, which is somewhat arbitrary, so create
          // a vector sorted by IDs for test purposes.
          std::vector<dof_id_type> neighbor_node_ids(neighbor_nodes.size());
          for (unsigned i=0; i<neighbor_nodes.size(); ++i)
            neighbor_node_ids[i] = neighbor_nodes[i]->id();
          std::sort(neighbor_node_ids.begin(), neighbor_node_ids.end());

          // Compare to validation_data
          for (unsigned j=0; j<neighbor_node_ids.size(); ++j)
            {
              CPPUNIT_ASSERT_EQUAL( validation_data[2*ctr + j], neighbor_node_ids[j] );
            }
        }
    }
  }

public:
  void setUp() {}

  void tearDown() {}

  void testEdge2()
  {
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

};


CPPUNIT_TEST_SUITE_REGISTRATION( NodalNeighborsTest );
