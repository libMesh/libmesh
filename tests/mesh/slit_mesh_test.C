// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>

#include "test_comm.h"

using namespace libMesh;

class SlitMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh with nodes overlapping
   * on opposite sides of an internal, "slit" edge is useable.
   */
public:
  CPPUNIT_TEST_SUITE( SlitMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

protected:

  Mesh* _mesh;

  void build_mesh()
  {
    _mesh = new Mesh(*TestCommWorld);

    /*
      (0,1)           (1,1)           (2,1)
        x---------------x---------------x
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        x---------------x---------------x
       (0,0)           (1,0)          (2,0)
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        x---------------x---------------x
       (0,-1)          (1,-1)         (2,-1)
     */

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0, 0.0), 0 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(0.0, 0.0), 6 );
    _mesh->add_point( Point(2.0, 0.0), 7 );
    _mesh->add_point( Point(2.0, 1.0), 8 );
    _mesh->add_point( Point(2.0,-1.0), 9 );

    {
      Elem* elem_top_left = _mesh->add_elem( new Quad4 );
      elem_top_left->set_node(0) = _mesh->node_ptr(0);
      elem_top_left->set_node(1) = _mesh->node_ptr(1);
      elem_top_left->set_node(2) = _mesh->node_ptr(2);
      elem_top_left->set_node(3) = _mesh->node_ptr(3);

      Elem* elem_bottom_left = _mesh->add_elem( new Quad4 );
      elem_bottom_left->set_node(0) = _mesh->node_ptr(4);
      elem_bottom_left->set_node(1) = _mesh->node_ptr(5);
      elem_bottom_left->set_node(2) = _mesh->node_ptr(6);
      elem_bottom_left->set_node(3) = _mesh->node_ptr(0);

      Elem* elem_top_right = _mesh->add_elem( new Quad4 );
      elem_top_right->set_node(0) = _mesh->node_ptr(1);
      elem_top_right->set_node(1) = _mesh->node_ptr(7);
      elem_top_right->set_node(2) = _mesh->node_ptr(8);
      elem_top_right->set_node(3) = _mesh->node_ptr(2);

      Elem* elem_bottom_right = _mesh->add_elem( new Quad4 );
      elem_bottom_right->set_node(0) = _mesh->node_ptr(5);
      elem_bottom_right->set_node(1) = _mesh->node_ptr(9);
      elem_bottom_right->set_node(2) = _mesh->node_ptr(7);
      elem_bottom_right->set_node(3) = _mesh->node_ptr(6);
    }

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    _mesh->prepare_for_use(true /*skip_renumber*/);
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
   }

  void testMesh()
  {
    // There'd better be 4 elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)4, _mesh->n_elem() );

    // There'd better still be a full 10 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)10, _mesh->n_nodes() );

    /* The middle nodes should still be distinct between the top and
     * bottom elements */
    if (_mesh->query_elem(0) && _mesh->query_elem(1))
      CPPUNIT_ASSERT( _mesh->elem(0)->node(1) != _mesh->elem(1)->node(2) );
    if (_mesh->query_elem(2) && _mesh->query_elem(3))
      CPPUNIT_ASSERT( _mesh->elem(2)->node(0) != _mesh->elem(3)->node(3) );

    /* The middle nodes should still be shared between left and right
     * elements on top and bottom */
    if (_mesh->query_elem(0) && _mesh->query_elem(2))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem(0)->node(1), _mesh->elem(2)->node(0) );
    if (_mesh->query_elem(1) && _mesh->query_elem(3))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem(1)->node(2), _mesh->elem(3)->node(3) );
  }

};

class SlitMeshRefinedMeshTest : public SlitMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  CPPUNIT_TEST_SUITE( SlitMeshRefinedMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
    this->build_mesh();

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
#endif
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 20 total and 16 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)20, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)16, _mesh->n_active_elem() );

    // We should have 28 nodes, not 25 or 26
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)28, _mesh->n_nodes() );
#endif
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshRefinedMeshTest );
