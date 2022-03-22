#include <libmesh/replicated_mesh.h>
#include <libmesh/point.h>
#include <libmesh/elem.h>
#include <libmesh/face_quad4.h>

#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

class GetBoundaryPointsTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that ReplicatedMesh::get_boundary_points returns
   * a vector of boundary points enclosing the domain
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( GetBoundaryPointsTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  ReplicatedMesh* _mesh;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld, 2);

    // (0,1)           (1,1)           (2,1) (3,1)       (4,1)       (5,1)       (6,1)
    // x---------------x---------------x     x-----------x-----------x-----------x
    // |               |               |     |           |           |           |
    // |               |               |     |           |           |           |
    // |               |               |     x-----------x-----------x-----------x
    // |               |               |     |(3,0.5)    |  void     |           |(6,0.5)
    // |               |               |     |           |           |           |
    // x---------------x---------------x     x-----------x-----------x-----------x
    // (0,0)           (1,0)          (2,0)  | (3,0)     |           |           |(6,0)
    // |               |               |     |           |           |           |
    // |               |               |     |           |           |           |
    // |               |               |     |           |           |           |
    // |               |               |     |           |           |           |
    // x---------------x---------------x     x-----------x-----------x-----------x
    // (0,-1)          (1,-1)         (2,-1) (3,-1)                              (6,-1)

    _mesh->add_point( Point(0.0, 0.0), 0 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(2.0, 0.0), 6 );
    _mesh->add_point( Point(2.0, 1.0), 7 );
    _mesh->add_point( Point(2.0,-1.0), 8 );

    _mesh->add_point( Point(3.0,-1.0),  9 );
    _mesh->add_point( Point(3.0, 0.0), 10 );
    _mesh->add_point( Point(3.0, 0.5), 11 );
    _mesh->add_point( Point(3.0, 1.0), 12 );
    _mesh->add_point( Point(4.0,-1.0), 13 );
    _mesh->add_point( Point(4.0, 0.0), 14 );
    _mesh->add_point( Point(4.0, 0.5), 15 );
    _mesh->add_point( Point(4.0, 1.0), 16 );
    _mesh->add_point( Point(5.0,-1.0), 17 );
    _mesh->add_point( Point(5.0, 0.0), 18 );
    _mesh->add_point( Point(5.0, 0.5), 19 );
    _mesh->add_point( Point(5.0, 1.0), 20 );
    _mesh->add_point( Point(6.0,-1.0), 21 );
    _mesh->add_point( Point(6.0, 0.0), 22 );
    _mesh->add_point( Point(6.0, 0.5), 23 );
    _mesh->add_point( Point(6.0, 1.0), 24 );

    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 0));
      elem->set_node(0) = _mesh->node_ptr(0);
      elem->set_node(1) = _mesh->node_ptr(1);
      elem->set_node(2) = _mesh->node_ptr(2);
      elem->set_node(3) = _mesh->node_ptr(3);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 1));
      elem->set_node(0) = _mesh->node_ptr(4);
      elem->set_node(1) = _mesh->node_ptr(5);
      elem->set_node(2) = _mesh->node_ptr(1);
      elem->set_node(3) = _mesh->node_ptr(0);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 2));
      elem->set_node(0) = _mesh->node_ptr(7);
      elem->set_node(1) = _mesh->node_ptr(2);
      elem->set_node(2) = _mesh->node_ptr(1);
      elem->set_node(3) = _mesh->node_ptr(6);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 3));
      elem->set_node(0) = _mesh->node_ptr(5);
      elem->set_node(1) = _mesh->node_ptr(8);
      elem->set_node(2) = _mesh->node_ptr(6);
      elem->set_node(3) = _mesh->node_ptr(1);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 4));
      elem->set_node(0) = _mesh->node_ptr(9);
      elem->set_node(1) = _mesh->node_ptr(13);
      elem->set_node(2) = _mesh->node_ptr(14);
      elem->set_node(3) = _mesh->node_ptr(10);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 5));
      elem->set_node(0) = _mesh->node_ptr(10);
      elem->set_node(1) = _mesh->node_ptr(14);
      elem->set_node(2) = _mesh->node_ptr(15);
      elem->set_node(3) = _mesh->node_ptr(11);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 6));
      elem->set_node(0) = _mesh->node_ptr(11);
      elem->set_node(1) = _mesh->node_ptr(15);
      elem->set_node(2) = _mesh->node_ptr(16);
      elem->set_node(3) = _mesh->node_ptr(12);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 7));
      elem->set_node(0) = _mesh->node_ptr(13);
      elem->set_node(1) = _mesh->node_ptr(17);
      elem->set_node(2) = _mesh->node_ptr(18);
      elem->set_node(3) = _mesh->node_ptr(14);
    }
    // skip one element here
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 8));
      elem->set_node(0) = _mesh->node_ptr(15);
      elem->set_node(1) = _mesh->node_ptr(19);
      elem->set_node(2) = _mesh->node_ptr(20);
      elem->set_node(3) = _mesh->node_ptr(16);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 9));
      elem->set_node(0) = _mesh->node_ptr(17);
      elem->set_node(1) = _mesh->node_ptr(21);
      elem->set_node(2) = _mesh->node_ptr(22);
      elem->set_node(3) = _mesh->node_ptr(18);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 10));
      elem->set_node(0) = _mesh->node_ptr(18);
      elem->set_node(1) = _mesh->node_ptr(22);
      elem->set_node(2) = _mesh->node_ptr(23);
      elem->set_node(3) = _mesh->node_ptr(19);
    }
    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 11));
      elem->set_node(0) = _mesh->node_ptr(19);
      elem->set_node(1) = _mesh->node_ptr(23);
      elem->set_node(2) = _mesh->node_ptr(24);
      elem->set_node(3) = _mesh->node_ptr(20);
    }

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    _mesh->allow_renumbering(false);

    _mesh->prepare_for_use();
  }

public:
  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();
#endif
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    const auto points = _mesh->get_boundary_points();

    // There'd better be two disconnected subdomains
    CPPUNIT_ASSERT_EQUAL( (std::size_t)2, points.size() );

    // The first key should better be element 0
    auto it = points.find(0);
    CPPUNIT_ASSERT(it != points.end());

    // There'd better be one boundary in the first subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)1, it->second.size() );

    // There'd better be eight points on the boundary of the first subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)8, it->second[0].size() );

    // Boundary points should better start from (1,1)
    CPPUNIT_ASSERT_EQUAL( Point(1,1), it->second[0][0] );

    // Boundary points should better end with (2,1)
    CPPUNIT_ASSERT_EQUAL( Point(2,1), it->second[0][7] );

    // The second key should better be element 4
    it = points.find(4);
    CPPUNIT_ASSERT(it != points.end());

    // There'd better be two boundaries in the second subdomain due to the middle hole
    CPPUNIT_ASSERT_EQUAL( (std::size_t)2, it->second.size() );

    // There'd better be 12 points on the first (outer) boundary of the second subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)12, it->second[0].size() );

    // There'd better be 4 points on the first (inner) boundary of the second subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)4, it->second[1].size() );

    // The first boundary points should better start from (3,-1)
    CPPUNIT_ASSERT_EQUAL( Point(3,-1), it->second[0][0] );

    // The first boundary points should better end with (3,0)
    CPPUNIT_ASSERT_EQUAL( Point(3,0), it->second[0][11] );

    // The second boundary points should better start from (4,0)
    CPPUNIT_ASSERT_EQUAL( Point(4,0), it->second[1][0] );

    // The second boundary points should better end with (5,0)
    CPPUNIT_ASSERT_EQUAL( Point(5,0), it->second[1][3] );
  }

};

class GetBoundaryPointsSecondTest : public GetBoundaryPointsTest {
  /**
   * The goal of this test is the same as the previous, but now we
   * use the second order mesh.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( GetBoundaryPointsSecondTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
#endif

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();
    _mesh->all_second_order(true);
#endif
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    const auto points = _mesh->get_boundary_points();

    // There'd better be two disconnected subdomains
    CPPUNIT_ASSERT_EQUAL( (std::size_t)2, points.size() );

    // The first key should better be element 0
    auto it = points.find(0);
    CPPUNIT_ASSERT(it != points.end());

    // There'd better be one boundary in the first subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)1, it->second.size() );

    // There'd better be 16 points on the boundary of the first subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)16, it->second[0].size() );

    // Boundary points should better start from (1,1)
    CPPUNIT_ASSERT_EQUAL( Point(1,1), it->second[0][0] );

    // Boundary points should better include side middle points including  (0.5,1)
    CPPUNIT_ASSERT_EQUAL( Point(0.5,1), it->second[0][1] );

    // Boundary points should better end with (1.5,1)
    CPPUNIT_ASSERT_EQUAL( Point(1.5,1), it->second[0][15] );

    // The second key should better be element 4
    it = points.find(4);
    CPPUNIT_ASSERT(it != points.end());

    // There'd better be two boundaries in the second subdomain due to the middle hole
    CPPUNIT_ASSERT_EQUAL( (std::size_t)2, it->second.size() );

    // There'd better be 24 points on the first (outer) boundary of the second subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)24, it->second[0].size() );

    // There'd better be 8 points on the first (inner) boundary of the second subdomain
    CPPUNIT_ASSERT_EQUAL( (std::size_t)8, it->second[1].size() );

    // The first boundary points should better start from (3,-1)
    CPPUNIT_ASSERT_EQUAL( Point(3,-1), it->second[0][0] );

    // The first boundary points should better include side middle points including (3.5,-1)
    CPPUNIT_ASSERT_EQUAL( Point(3.5,-1), it->second[0][1] );

    // The first boundary points should better end with (3,-0.5)
    CPPUNIT_ASSERT_EQUAL( Point(3,-0.5), it->second[0][23] );

    // The second boundary points should better start from (4,0)
    CPPUNIT_ASSERT_EQUAL( Point(4,0), it->second[1][0] );

    // The second boundary points should better include side middle points including (4,0.25)
    CPPUNIT_ASSERT_EQUAL( Point(4,0.25), it->second[1][1] );

    // The second boundary points should better end with (4.5,0)
    CPPUNIT_ASSERT_EQUAL( Point(4.5,0), it->second[1][7] );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( GetBoundaryPointsTest );
CPPUNIT_TEST_SUITE_REGISTRATION( GetBoundaryPointsSecondTest );
