// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/serial_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>

#include "test_comm.h"

using namespace libMesh;

class MixedDimensionMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh with 1D elements overlapping
   * on the edge is consistent. That is, they share the same global node numbers and
   * the same dof numbers for a variable.
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionMeshTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );
  CPPUNIT_TEST( testPointLocatorList );
  CPPUNIT_TEST( testPointLocatorTree );

  CPPUNIT_TEST_SUITE_END();

protected:

  SerialMesh* _mesh;

  void build_mesh()
  {
    _mesh = new SerialMesh(*TestCommWorld);

    /*
      (0,1)           (1,1)
        x---------------x
        |               |
        |               |
        |               |
        |               |
        |               |
        x---------------x
       (0,0)           (1,0)
        |               |
        |               |
        |               |
        |               |
        x---------------x
       (0,-1)          (1,-1)
     */

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0, 0.0), 0 );

    {
      Elem* elem_top = _mesh->add_elem( new Quad4 );
      elem_top->set_node(0) = _mesh->node_ptr(0);
      elem_top->set_node(1) = _mesh->node_ptr(1);
      elem_top->set_node(2) = _mesh->node_ptr(2);
      elem_top->set_node(3) = _mesh->node_ptr(3);

      Elem* elem_bottom = _mesh->add_elem( new Quad4 );
      elem_bottom->set_node(0) = _mesh->node_ptr(4);
      elem_bottom->set_node(1) = _mesh->node_ptr(5);
      elem_bottom->set_node(2) = _mesh->node_ptr(1);
      elem_bottom->set_node(3) = _mesh->node_ptr(0);

      Elem* edge = _mesh->add_elem( new Edge2 );
      edge->set_node(0) = _mesh->node_ptr(0);
      edge->set_node(1) = _mesh->node_ptr(1);

      // 2D elements will have subdomain id 0, this one will have 1
      edge->subdomain_id() = 1;
    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequenly,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);
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
    // There'd better be 3 elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)3, _mesh->n_elem() );

    // There'd better be only 6 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)6, _mesh->n_nodes() );

    /* The nodes for the EDGE2 element should have the same global ids
       as the bottom edge of the top QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(2)->node(0), _mesh->elem(0)->node(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(2)->node(1), _mesh->elem(0)->node(1) );

    /* The nodes for the EDGE2 element should have the same global ids
       as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(2)->node(0), _mesh->elem(1)->node(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(2)->node(1), _mesh->elem(1)->node(2) );

    /* The nodes for the bottom edge of the top QUAD4 element should have
       the same global ids as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(0)->node(0), _mesh->elem(1)->node(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(0)->node(1), _mesh->elem(1)->node(2) );

    // We didn't set an interior_parent on the edge element, so it
    // should default to NULL
    CPPUNIT_ASSERT( !_mesh->elem(2)->interior_parent() );
  }

  void testDofOrdering()
  {
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    std::vector<dof_id_type> top_quad_dof_indices, bottom_quad_dof_indices, edge_dof_indices;

    dof_map.dof_indices( _mesh->elem(0), top_quad_dof_indices );
    dof_map.dof_indices( _mesh->elem(1), bottom_quad_dof_indices );
    dof_map.dof_indices( _mesh->elem(2), edge_dof_indices );

    /* The dofs for the EDGE2 element should be the same
       as the bottom edge of the top QUAD4 dofs */
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[0], top_quad_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[1], top_quad_dof_indices[1] );

    /* The dofs for the EDGE2 element should be the same
       as the top edge of the bottom QUAD4 dofs */
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[0], bottom_quad_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[1], bottom_quad_dof_indices[2] );

    /* The nodes for the bottom edge of the top QUAD4 element should have
       the same global ids as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( top_quad_dof_indices[0], bottom_quad_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( top_quad_dof_indices[1], bottom_quad_dof_indices[2] );
  }

  void testPointLocatorList()
  {
    UniquePtr<PointLocatorBase> locator = PointLocatorBase::build(LIST,*_mesh);

    Point top_point(0.4, 0.5);
    const Elem* top_elem = (*locator)(top_point);
    CPPUNIT_ASSERT(top_elem);

    // We should have gotten back the top quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)0, top_elem->id() );

    Point bottom_point(0.5, -0.5);
    const Elem* bottom_elem = (*locator)(bottom_point);
    CPPUNIT_ASSERT(bottom_elem);

    // We should have gotten back the bottom quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)1, bottom_elem->id() );

    // Test getting back the edge
    {
      std::set<subdomain_id_type> subdomain_id; subdomain_id.insert(1);
      Point interface_point( 0.2, 0.0 );
      const Elem* interface_elem = (*locator)(interface_point, &subdomain_id);
      CPPUNIT_ASSERT(interface_elem);

      // We should have gotten back the overlapping edge element
      CPPUNIT_ASSERT_EQUAL( (dof_id_type)2, interface_elem->id() );
    }
  }

  void testPointLocatorTree()
  {
    UniquePtr<PointLocatorBase> locator = _mesh->sub_point_locator();

    Point top_point(0.5, 0.5);
    const Elem* top_elem = (*locator)(top_point);
    CPPUNIT_ASSERT(top_elem);

    // We should have gotten back the top quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)0, top_elem->id() );

    Point bottom_point(0.5, -0.5);
    const Elem* bottom_elem = (*locator)(bottom_point);
    CPPUNIT_ASSERT(bottom_elem);

    // We should have gotten back the bottom quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)1, bottom_elem->id() );

    // Test getting back the edge
    {
      std::set<subdomain_id_type> subdomain_id; subdomain_id.insert(1);
      Point interface_point( 0.5, 0.0 );
      const Elem* interface_elem = (*locator)(interface_point, &subdomain_id);
      CPPUNIT_ASSERT(interface_elem);

      // We should have gotten back the overlapping edge element
      CPPUNIT_ASSERT_EQUAL( (dof_id_type)2, interface_elem->id() );
    }
  }

};

class MixedDimensionRefinedMeshTest : public MixedDimensionMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionRefinedMeshTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
    /*

        3-------10------2
        |       |       |
        |   5   |   6   |
        8-------7-------9
        |       |       |
        |   3   |   4   |
        0-------6-------1
        |       |       |
        |   9   |  10   |
       13------12-------14
        |       |       |
        |   7   |   8   |
        4-------11------5

     */
    this->build_mesh();

    // Let's set an interior_parent() this time for testing
    _mesh->elem(2)->set_interior_parent(_mesh->elem(0));

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
#endif
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 13 total and 10 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)13, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)10, _mesh->n_active_elem() );

    // We should have 15 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)15, _mesh->n_nodes() );

    // EDGE2,id=11 should have same nodes of bottom of QUAD4, id=3
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->node(0), _mesh->elem(3)->node(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->node(1), _mesh->elem(3)->node(1) );

    // EDGE2,id=12 should have same nodes of bottom of QUAD4, id=4
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->node(0), _mesh->elem(4)->node(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->node(1), _mesh->elem(4)->node(1) );

    // EDGE2,id=11 should have same nodes of top of QUAD4, id=9
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->node(0), _mesh->elem(9)->node(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->node(1), _mesh->elem(9)->node(2) );

    // EDGE2,id=12 should have same nodes of top of QUAD4, id=10
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->node(0), _mesh->elem(10)->node(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->node(1), _mesh->elem(10)->node(2) );

    // Shared node between the EDGE2 elements should have the same global id
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->node(1), _mesh->elem(12)->node(0) );

    // EDGE2 child elements should have the correct parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->parent(), _mesh->elem(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->parent(), _mesh->elem(2) );

    // EDGE2 child elements should have the correct interior_parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(11)->interior_parent(), _mesh->elem(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem(12)->interior_parent(), _mesh->elem(4) );
#endif
  }

  void testDofOrdering()
  {
#ifdef LIBMESH_ENABLE_AMR
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    std::vector<dof_id_type> top_quad3_dof_indices, top_quad4_dof_indices;
    std::vector<dof_id_type> bottom_quad9_dof_indices, bottom_quad10_dof_indices;
    std::vector<dof_id_type> edge11_dof_indices, edge12_dof_indices;

    dof_map.dof_indices( _mesh->elem(3), top_quad3_dof_indices );
    dof_map.dof_indices( _mesh->elem(4), top_quad4_dof_indices );
    dof_map.dof_indices( _mesh->elem(9), bottom_quad9_dof_indices );
    dof_map.dof_indices( _mesh->elem(10), bottom_quad10_dof_indices );
    dof_map.dof_indices( _mesh->elem(11), edge11_dof_indices );
    dof_map.dof_indices( _mesh->elem(12), edge12_dof_indices );

    // EDGE2,id=11 should have same dofs as of bottom of QUAD4, id=3
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[0], top_quad3_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], top_quad3_dof_indices[1] );

    // EDGE2,id=12 should have same dofs of bottom of QUAD4, id=4
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[0], top_quad4_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[1], top_quad4_dof_indices[1] );

    // EDGE2,id=11 should have same dofs of top of QUAD4, id=9
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[0], bottom_quad9_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], bottom_quad9_dof_indices[2] );

    // EDGE2,id=12 should have same dofs of top of QUAD4, id=10
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[0], bottom_quad10_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[1], bottom_quad10_dof_indices[2] );

    //EDGE2 elements should have same shared dof number
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], edge12_dof_indices[0] );
#endif
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionRefinedMeshTest );
