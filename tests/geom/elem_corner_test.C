#include <libmesh/elem_corner.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

#if LIBMESH_DIM > 1
class ElemCornerTest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( ElemCornerTest );
  CPPUNIT_TEST( test );
  CPPUNIT_TEST_SUITE_END();

public:

  void test()
  {
    LOG_UNIT_TEST;

    // QUAD elem for vertex testing
    Mesh mesh_2d (*TestCommWorld);
    MeshTools::Generation::build_square (mesh_2d, 1, 1);
    const auto & quad = mesh_2d.elem_ref(0);

    // Default constructor: unvalid state
    const ElemCorner default_corner;
    CPPUNIT_ASSERT(default_corner.is_invalid());
    CPPUNIT_ASSERT(!default_corner.at_corner());
    CPPUNIT_ASSERT(!default_corner.at_vertex());
#if LIBMESH_DIM > 2
    CPPUNIT_ASSERT(!default_corner.at_edge());
#endif // LIBMESH_DIM > 2
    CPPUNIT_ASSERT_EQUAL((std::string)"not at corner", default_corner.print());
    CPPUNIT_ASSERT(default_corner.is_valid(quad, quad.point(0)));

    // At vertex
    const ElemCorner at_vertex(0, Elem::invalid_vertex);
    CPPUNIT_ASSERT(!at_vertex.is_invalid());
    CPPUNIT_ASSERT(at_vertex.at_corner());
    CPPUNIT_ASSERT(at_vertex.at_vertex());
#if LIBMESH_DIM> 2
    CPPUNIT_ASSERT(!at_vertex.at_edge());
    CPPUNIT_ASSERT(!at_vertex.at_edge(0, 1));
    CPPUNIT_ASSERT(!at_vertex.at_edge(quad, 0));
#endif // LIBMESH_DIM > 2
    CPPUNIT_ASSERT(at_vertex.at_vertex(0));
    CPPUNIT_ASSERT(!at_vertex.at_vertex(1));
    CPPUNIT_ASSERT_EQUAL(at_vertex.vertex(), (unsigned short)0);
    CPPUNIT_ASSERT(at_vertex.is_valid(quad, quad.point(0)));
    CPPUNIT_ASSERT(!at_vertex.is_valid(quad, quad.point(2)));

#if LIBMESH_DIM > 2
    // HEX elem for edge testing
    Mesh mesh_3d (*TestCommWorld);
    MeshTools::Generation::build_cube (mesh_3d, 1, 1, 1);
    const auto & hex = mesh_3d.elem_ref(0);

    // At edge
    const ElemCorner at_edge(0, 1);
    CPPUNIT_ASSERT(!at_edge.is_invalid());
    CPPUNIT_ASSERT(at_edge.at_corner());
    CPPUNIT_ASSERT(!at_edge.at_vertex());
    CPPUNIT_ASSERT(!at_edge.at_vertex(0));
    CPPUNIT_ASSERT(at_edge.at_edge());
    CPPUNIT_ASSERT(at_edge.at_edge(0, 1));
    CPPUNIT_ASSERT(at_edge.at_edge(1, 0));
    CPPUNIT_ASSERT(at_edge.at_edge(hex, 0));
    CPPUNIT_ASSERT_EQUAL(at_edge.edge_vertices().first, (unsigned short)0);
    CPPUNIT_ASSERT_EQUAL(at_edge.edge_vertices().second, (unsigned short)1);
    CPPUNIT_ASSERT(!at_edge.is_invalid());
    CPPUNIT_ASSERT(at_edge.is_valid(hex, (hex.point(0) + hex.point(1)) * 0.5));
    CPPUNIT_ASSERT(at_edge.is_valid(hex, hex.point(0)));
    CPPUNIT_ASSERT(at_edge.is_valid(hex, hex.point(1)));
    CPPUNIT_ASSERT(!at_edge.is_valid(hex, (hex.point(2)) * 0.5));

    // build_edge
    auto edge_elem = at_edge.build_edge(hex);
    CPPUNIT_ASSERT((edge_elem->point(0).absolute_fuzzy_equals(hex.point(0))
                    && edge_elem->point(1).absolute_fuzzy_equals(hex.point(1))) ||
                   ((edge_elem->point(1).absolute_fuzzy_equals(hex.point(0))
                    && edge_elem->point(0).absolute_fuzzy_equals(hex.point(1)))));
#endif // LIBMESH_DIM > 2

    // For testing non-const methods
    ElemCorner corner(0, 0);

    // invalidate
    corner.invalidate();
    CPPUNIT_ASSERT(corner.is_invalid());

    // set_vertex
    corner.set_vertex(0);
    CPPUNIT_ASSERT(corner.at_vertex(0));

    // set_edge
#if LIBMESH_DIM > 2
    corner.set_edge(0, 1);
    CPPUNIT_ASSERT(corner.at_edge(0, 1));
    corner.set_edge(std::pair<unsigned short, unsigned short>(1, 2));
    CPPUNIT_ASSERT(corner.at_edge(1, 2));
#endif // LIBMESH_DIM > 2
  }

};
#endif // LIBMESH_DIM > 1

CPPUNIT_TEST_SUITE_REGISTRATION( ElemCornerTest );
