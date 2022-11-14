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

    // Default constructor: unvalid state
    const ElemCorner default_corner;
    CPPUNIT_ASSERT(default_corner.is_invalid());
    CPPUNIT_ASSERT(!default_corner.at_corner());
    CPPUNIT_ASSERT(!default_corner.at_vertex());
    CPPUNIT_ASSERT(!default_corner.at_edge());
    CPPUNIT_ASSERT_EQUAL((std::string)"not at corner", default_corner.print());

    // EDGE elem for vertex testing
    Mesh mesh_1d (*TestCommWorld);
    MeshTools::Generation::build_line(mesh_1d, 1);
    const auto edge = mesh_1d.elem_ptr(0);

    // At edge vertex
    const ElemCorner at_edge_vertex(1, Elem::invalid_vertex);
    CPPUNIT_ASSERT(!at_edge_vertex.is_invalid());
    CPPUNIT_ASSERT(at_edge_vertex.at_corner());
    CPPUNIT_ASSERT(at_edge_vertex.at_vertex());
    CPPUNIT_ASSERT(!at_edge_vertex.at_edge());
    CPPUNIT_ASSERT(at_edge_vertex.at_vertex(1));
    CPPUNIT_ASSERT_EQUAL(at_edge_vertex.vertex(), (unsigned short)1);
    if (edge)
      for (const auto v : edge->vertex_index_range())
      {
        CPPUNIT_ASSERT_EQUAL(at_edge_vertex.at_vertex(v), v == 1);
        CPPUNIT_ASSERT_EQUAL(at_edge_vertex.is_valid(*edge, edge->point(v)), v == 1);
      }

#if LIBMESH_DIM > 1
    // QUAD elem for vertex testing
    Mesh mesh_2d (*TestCommWorld);
    MeshTools::Generation::build_square (mesh_2d, 1, 1);
    const auto quad = mesh_2d.elem_ptr(0);

    // At quad vertex
    const ElemCorner at_quad_vertex(0, Elem::invalid_vertex);
    CPPUNIT_ASSERT(!at_quad_vertex.is_invalid());
    CPPUNIT_ASSERT(at_quad_vertex.at_corner());
    CPPUNIT_ASSERT(at_quad_vertex.at_vertex());
    CPPUNIT_ASSERT(!at_quad_vertex.at_edge());
    CPPUNIT_ASSERT(!at_quad_vertex.at_edge(0, 1));
    CPPUNIT_ASSERT(at_quad_vertex.at_vertex(0));
    CPPUNIT_ASSERT(!at_quad_vertex.at_vertex(1));
    CPPUNIT_ASSERT_EQUAL(at_quad_vertex.vertex(), (unsigned short)0);
    if (quad)
    {
      for (const auto v : quad->vertex_index_range())
        CPPUNIT_ASSERT(default_corner.is_valid(*quad, quad->point(v)));
      for (const auto e : quad->edge_index_range())
        CPPUNIT_ASSERT(!at_quad_vertex.at_edge(*quad, e));
      for (const auto v : quad->vertex_index_range())
        CPPUNIT_ASSERT_EQUAL(at_quad_vertex.is_valid(*quad, quad->point(v)), v == 0);
    }

    // At quad "edge" (a side)
    const ElemCorner at_quad_edge(1, 2);
    CPPUNIT_ASSERT(!at_quad_edge.is_invalid());
    CPPUNIT_ASSERT(at_quad_edge.at_corner());
    CPPUNIT_ASSERT(!at_quad_edge.at_vertex());
    CPPUNIT_ASSERT(at_quad_edge.at_edge());
    CPPUNIT_ASSERT(!at_quad_edge.at_edge(0, 1));
    CPPUNIT_ASSERT(at_quad_edge.at_edge(1, 2));
    CPPUNIT_ASSERT_EQUAL(at_quad_edge.edge_vertices().first, (unsigned short)1);
    CPPUNIT_ASSERT_EQUAL(at_quad_edge.edge_vertices().second, (unsigned short)2);
    if (quad)
    {
      CPPUNIT_ASSERT(at_quad_edge.is_valid(*quad, 0.5 * (quad->point(1) + quad->point(2))));
      for (const auto e : quad->edge_index_range())
        CPPUNIT_ASSERT_EQUAL(at_quad_edge.at_edge(*quad, e), e == 1);
    }
#endif

#if LIBMESH_DIM > 2
    // HEX elem for edge testing
    Mesh mesh_3d (*TestCommWorld);
    MeshTools::Generation::build_cube (mesh_3d, 1, 1, 1);
    const auto hex = mesh_3d.elem_ptr(0);

    // At edge
    const ElemCorner at_hex_edge(0, 1);
    CPPUNIT_ASSERT(!at_hex_edge.is_invalid());
    CPPUNIT_ASSERT(at_hex_edge.at_corner());
    CPPUNIT_ASSERT(!at_hex_edge.at_vertex());
    CPPUNIT_ASSERT(!at_hex_edge.at_vertex(0));
    CPPUNIT_ASSERT(at_hex_edge.at_edge());
    CPPUNIT_ASSERT(at_hex_edge.at_edge(0, 1));
    CPPUNIT_ASSERT(at_hex_edge.at_edge(1, 0));
    CPPUNIT_ASSERT_EQUAL(at_hex_edge.edge_vertices().first, (unsigned short)0);
    CPPUNIT_ASSERT_EQUAL(at_hex_edge.edge_vertices().second, (unsigned short)1);
    CPPUNIT_ASSERT(!at_hex_edge.is_invalid());
    if (hex)
    {
      CPPUNIT_ASSERT(at_hex_edge.is_valid(*hex, (hex->point(0) + hex->point(1)) * 0.5));
      for (const auto e : hex->edge_index_range())
        CPPUNIT_ASSERT_EQUAL(at_hex_edge.at_edge(*hex, e), e == 0);
      for (const auto v : hex->vertex_index_range())
        CPPUNIT_ASSERT_EQUAL(at_hex_edge.is_valid(*hex, hex->point(v)), v == 0 || v == 1);
      CPPUNIT_ASSERT(!at_hex_edge.is_valid(*hex, (hex->point(2)) * 0.5));

      // build_edge
      auto edge_elem = at_hex_edge.build_edge(*hex);
      CPPUNIT_ASSERT((edge_elem->point(0).absolute_fuzzy_equals(hex->point(0))
                      && edge_elem->point(1).absolute_fuzzy_equals(hex->point(1))) ||
                  ((edge_elem->point(1).absolute_fuzzy_equals(hex->point(0))
                       && edge_elem->point(0).absolute_fuzzy_equals(hex->point(1)))));
    }
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
