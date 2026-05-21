#include <libmesh/libmesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/cell_c0polyhedron.h>
#include <libmesh/cell_polyhedron.h>
#include <libmesh/face_c0polygon.h>
#include <libmesh/face_polygon.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>


using namespace libMesh;

class AllTriTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the Mesh Extruder
   * with the optional object callback for setting custom subdomain IDs.
   * We pass a custom object for generating subdomains based on the old element
   * ID and the current layer and assert the proper values.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllTriTest );

  // 2D tests
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testAllTriTri );
  CPPUNIT_TEST( testAllTriQuad );
  CPPUNIT_TEST( testAllTriQuad8 );
  CPPUNIT_TEST( testAllTriQuad9 );
  CPPUNIT_TEST( testAllTriC0Polygon );
  CPPUNIT_TEST( testAllTriC0PolygonOctagon );
#endif

  // 3D tests
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testAllTriPrism6 );
  CPPUNIT_TEST( testAllTriPrism18 );
  CPPUNIT_TEST( testAllTriPrism20 );
  CPPUNIT_TEST( testAllTriPrism21 );
  CPPUNIT_TEST( testAllTriC0PolyhedronCube );
  CPPUNIT_TEST( testAllTriC0PolyhedronHexagonalPrism );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  // Helper function called by the test implementations, saves a few lines of code.
  void test_helper_2D(ElemType elem_type,
                      dof_id_type n_elem_expected,
                      std::size_t n_boundary_conds_expected)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    // Build a 2x1 TRI3 mesh and ask to split it into triangles.
    // Should be a no-op
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/2, /*ny=*/1,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        elem_type);

    MeshTools::Modification::all_tri(mesh);

    // Make sure that the expected number of elements is found.
    CPPUNIT_ASSERT_EQUAL(n_elem_expected, mesh.n_elem());

    // Make sure the expected number of BCs is found.
    CPPUNIT_ASSERT_EQUAL(n_boundary_conds_expected, mesh.get_boundary_info().n_boundary_conds());
  }

  // Helper function called by the test implementations in 3D, saves a few lines of code.
  void test_helper_3D(ElemType elem_type,
                      dof_id_type n_elem_expected,
                      std::size_t n_boundary_conds_expected)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/3);

    // Build a 2x1 TRI3 mesh and ask to split it into triangles.
    // Should be a no-op
    MeshTools::Generation::build_cube(mesh,
                                      /*nx=*/1, /*ny=*/1, /*nz=*/1,
                                      /*xmin=*/0., /*xmax=*/1.,
                                      /*ymin=*/0., /*ymax=*/1.,
                                      /*zmin=*/0., /*zmax=*/1.,
                                      elem_type);

    MeshTools::Modification::all_tri(mesh);

    // Make sure that the expected number of elements is found.
    CPPUNIT_ASSERT_EQUAL(n_elem_expected, mesh.n_elem());

    // Make sure the expected number of BCs is found.
    CPPUNIT_ASSERT_EQUAL(n_boundary_conds_expected, mesh.get_boundary_info().n_boundary_conds());
  }

public:
  void setUp() {}

  void tearDown() {}

  // 4 TRIs no-op
  void testAllTriTri() { LOG_UNIT_TEST; test_helper_2D(TRI3, /*nelem=*/4, /*nbcs=*/6); }

  // 2 quads split into 4 TRIs.
  void testAllTriQuad() { LOG_UNIT_TEST; test_helper_2D(QUAD4, /*nelem=*/4, /*nbcs=*/6); }

  // 2 QUAD8s split into 4 TRIs.
  void testAllTriQuad8() { LOG_UNIT_TEST; test_helper_2D(QUAD8, /*nelem=*/4, /*nbcs=*/6); }

  // 2 QUAD9s split into 4 TRIs.
  void testAllTriQuad9() { LOG_UNIT_TEST; test_helper_2D(QUAD9, /*nelem=*/4, /*nbcs=*/6); }

  // 2 PRISMs split into 6 TETs with 2 boundary faces per side.
  void testAllTriPrism6() { LOG_UNIT_TEST; test_helper_3D(PRISM6, /*nelem=*/6, /*nbcs=*/12); }
  void testAllTriPrism18() { LOG_UNIT_TEST; test_helper_3D(PRISM18, /*nelem=*/6, /*nbcs=*/12); }
  void testAllTriPrism20() { LOG_UNIT_TEST; test_helper_3D(PRISM20, /*nelem=*/6, /*nbcs=*/12); }
  void testAllTriPrism21() { LOG_UNIT_TEST; test_helper_3D(PRISM21, /*nelem=*/6, /*nbcs=*/12); }

  // Build a C0Polygon paving (triangles, quads, hexagons) via
  // build_square and split it into a pure TRI3 mesh.
  void testAllTriC0Polygon()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/2, /*ny=*/2,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        C0POLYGON);

    // The post-all_tri element count is the sum of the per-polygon
    // subtriangle counts, and external sides should be preserved one
    // for one as TRI3 sides.
    dof_id_type n_elem_expected = 0;
    for (const Elem * elem : mesh.element_ptr_range())
      {
        const Polygon * poly = dynamic_cast<const Polygon *>(elem);
        CPPUNIT_ASSERT(poly != nullptr);
        n_elem_expected += poly->n_subtriangles();
      }

    const std::size_t n_bcs_before =
      mesh.get_boundary_info().n_boundary_conds();

    MeshTools::Modification::all_tri(mesh);

    CPPUNIT_ASSERT_EQUAL(n_elem_expected, mesh.n_elem());
    CPPUNIT_ASSERT_EQUAL(n_bcs_before,
                         mesh.get_boundary_info().n_boundary_conds());

    for (const Elem * elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(ElemType(TRI3), elem->type());
  }

  // A single regular octagon (8 sides, 6 subtriangles) exercises the
  // path where a polygon requires more subelements than any non-polygon
  // 2D element type would.
  void testAllTriC0PolygonOctagon()
  {
    LOG_UNIT_TEST;

    constexpr unsigned int n_sides = 8;

    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    std::unique_ptr<Elem> octagon = std::make_unique<C0Polygon>(n_sides);
    for (unsigned int i = 0; i < n_sides; ++i)
      {
        const Real angle = 2 * libMesh::pi * i / n_sides;
        Node * node = mesh.add_point(Point(std::cos(angle), std::sin(angle), 0.),
                                     /*id=*/i);
        octagon->set_node(i, node);
      }
    octagon->set_id() = 0;
    Elem * elem = mesh.add_elem(std::move(octagon));

    // Mark every external side with a boundary id so we can verify
    // boundary information is transferred to the new triangles.
    for (unsigned int s = 0; s < n_sides; ++s)
      mesh.get_boundary_info().add_side(elem, s, /*bnd_id=*/0);

    mesh.prepare_for_use();

    MeshTools::Modification::all_tri(mesh);

    // n_sides - 2 = 6 subtriangles
    CPPUNIT_ASSERT_EQUAL(dof_id_type(n_sides - 2), mesh.n_elem());
    CPPUNIT_ASSERT_EQUAL(std::size_t(n_sides),
                         mesh.get_boundary_info().n_boundary_conds());

    for (const Elem * e : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(ElemType(TRI3), e->type());
  }

protected:

  // Builds a C0Polyhedron in mesh whose faces are described by
  // nodes_on_side (indices into the existing mesh node list), runs
  // all_tri, and verifies that the result is a pure TET4 mesh with the
  // expected sub-element count and preserved boundary data.
  void test_helper_c0polyhedron
    (const std::vector<Point> & points,
     const std::vector<std::vector<unsigned int>> & nodes_on_side)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/3);

    for (auto p : index_range(points))
      mesh.add_point(points[p], /*id=*/p);

    std::vector<std::shared_ptr<Polygon>> sides(nodes_on_side.size());
    for (auto s : index_range(nodes_on_side))
      {
        const auto & nodes_on_s = nodes_on_side[s];
        sides[s] = std::make_shared<C0Polygon>(nodes_on_s.size());
        for (auto i : index_range(nodes_on_s))
          sides[s]->set_node(i, mesh.node_ptr(nodes_on_s[i]));
      }

    std::unique_ptr<Node> mid_elem_node;
    std::unique_ptr<Elem> polyhedron =
      std::make_unique<C0Polyhedron>(sides, mid_elem_node);
    if (mid_elem_node)
      mesh.add_node(std::move(mid_elem_node));
    polyhedron->set_id() = 0;
    Elem * elem = mesh.add_elem(std::move(polyhedron));

    const auto * poly = cast_ptr<const C0Polyhedron *>(elem);
    const dof_id_type n_elem_expected = poly->n_subelements();

    // Mark every external face with a boundary id so we can verify
    // boundary information is transferred to the new tets.
    for (unsigned int s = 0; s < elem->n_sides(); ++s)
      mesh.get_boundary_info().add_side(elem, s, /*bnd_id=*/0);

    // The number of boundary triangles produced is the total number of
    // subtriangles across the polyhedron's polygonal faces.
    std::size_t n_bcs_expected = 0;
    for (unsigned int s = 0; s < elem->n_sides(); ++s)
      n_bcs_expected += sides[s]->n_subtriangles();

    mesh.prepare_for_use();

    MeshTools::Modification::all_tri(mesh);

    CPPUNIT_ASSERT_EQUAL(n_elem_expected, mesh.n_elem());
    CPPUNIT_ASSERT_EQUAL(n_bcs_expected,
                         mesh.get_boundary_info().n_boundary_conds());

    for (const Elem * e : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(ElemType(TET4), e->type());
  }

public:

  // A cube built as a C0Polyhedron exercises the path where the
  // optimal tetrahedralization succeeds (no mid-element node needed).
  void testAllTriC0PolyhedronCube()
  {
    LOG_UNIT_TEST;

    const std::vector<Point> points =
      { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };

    const std::vector<std::vector<unsigned int>> nodes_on_side =
      { {0, 1, 2, 3},   // min z
        {0, 1, 5, 4},   // min y
        {2, 6, 5, 1},   // max x
        {2, 3, 7, 6},   // max y
        {0, 4, 7, 3},   // min x
        {5, 6, 7, 4} }; // max z

    test_helper_c0polyhedron(points, nodes_on_side);
  }

  // A hexagonal prism exercises the fallback path where a mid-element
  // node is added to tetrahedralize the polyhedron.
  void testAllTriC0PolyhedronHexagonalPrism()
  {
    LOG_UNIT_TEST;

    const std::vector<Point> points =
      { { 0, -2, 0}, {-1, -1, 0}, {-1, 1, 0},
        { 0,  2, 0}, { 1,  1, 0}, { 1, -1, 0},
        { 0, -2, 1}, {-1, -1, 1}, {-1, 1, 1},
        { 0,  2, 1}, { 1,  1, 1}, { 1, -1, 1} };

    const std::vector<std::vector<unsigned int>> nodes_on_side =
      { {0, 1, 2, 3, 4, 5},
        {0, 1,  7,  6},
        {1, 2,  8,  7},
        {2, 3,  9,  8},
        {3, 4, 10,  9},
        {4, 5, 11, 10},
        {5, 0,  6, 11},
        {6, 7,  8,  9, 10, 11} };

    test_helper_c0polyhedron(points, nodes_on_side);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllTriTest );
