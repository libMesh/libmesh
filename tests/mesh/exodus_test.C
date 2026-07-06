#include "mesh_elem_test.h"

#ifdef LIBMESH_HAVE_EXODUS_API

#include "libmesh/cell_c0polyhedron.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/face_c0polygon.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/node.h"

#include <memory>
#include <string>
#include <vector>

using namespace libMesh;

template <ElemType elem_type>
class ExodusTest : public MeshPerElemTest<elem_type>
{
public:
  void test_read_gold()
  {
    LOG_UNIT_TEST;

    Mesh input_mesh(*TestCommWorld);

    ExodusII_IO exii(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii.read("meshes/exodus_elements/read_exodus_" +
                Utility::enum_to_string(elem_type) + ".e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT(this->meshes_equal_enough(input_mesh, true));
  }

  void test_write()
  {
    LOG_UNIT_TEST;

    // This is a *buffered* write; we use scope to make sure the
    // ExodusII_IO object gets destructed (and thus is guaranteed to
    // finish writing and close the file) before we try to read what
    // was written.
    {
      ExodusII_IO exii(*this->_mesh);

      // We still default to 32-char names for backwards
      // compatibility, but we're writing a mesh with extra-long names
      // in it for testing, so we manually enable longer names.
      exii.set_max_name_length(80);

      exii.write("write_exodus_" +
                 Utility::enum_to_string(elem_type) + ".e");
    }

    Mesh input_mesh(*TestCommWorld);
    ExodusII_IO exii_input(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii_input.read("write_exodus_" +
                      Utility::enum_to_string(elem_type) + ".e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT(this->meshes_equal_enough(input_mesh, true));
  }
};

#define EXODUSTEST              \
  CPPUNIT_TEST(test_read_gold); \
  CPPUNIT_TEST(test_write);

#define INSTANTIATE_EXODUSTEST(elemtype)                    \
  class ExodusTest_##elemtype : public ExodusTest<elemtype> \
  {                                                         \
  public:                                                   \
    ExodusTest_##elemtype() : ExodusTest<elemtype>()        \
    {                                                       \
      if (unitlog->summarized_logs_enabled())               \
        this->libmesh_suite_name = "ExodusTest";            \
      else                                                  \
        this->libmesh_suite_name = "ExodusTest_" #elemtype; \
    }                                                       \
    CPPUNIT_TEST_SUITE(ExodusTest_##elemtype);              \
    EXODUSTEST;                                             \
    CPPUNIT_TEST_SUITE_END();                               \
  };                                                        \
                                                            \
  CPPUNIT_TEST_SUITE_REGISTRATION(ExodusTest_##elemtype)

INSTANTIATE_EXODUSTEST(EDGE2);
INSTANTIATE_EXODUSTEST(EDGE3);
INSTANTIATE_EXODUSTEST(EDGE4);

#if LIBMESH_DIM > 1
INSTANTIATE_EXODUSTEST(TRI3);
INSTANTIATE_EXODUSTEST(TRISHELL3);
INSTANTIATE_EXODUSTEST(TRI6);
INSTANTIATE_EXODUSTEST(TRI7);

INSTANTIATE_EXODUSTEST(QUAD4);
INSTANTIATE_EXODUSTEST(QUADSHELL4);
INSTANTIATE_EXODUSTEST(QUAD8);
INSTANTIATE_EXODUSTEST(QUADSHELL8);
INSTANTIATE_EXODUSTEST(QUAD9);
INSTANTIATE_EXODUSTEST(QUADSHELL9);

class ExodusC0PolygonTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(ExodusC0PolygonTest);

  CPPUNIT_TEST(test_write_and_read_pentagon);

  CPPUNIT_TEST_SUITE_END();

  void build_pentagon(Mesh &mesh)
  {
    const std::vector<Point> points =
      { {0, 0}, {1, 0}, {1.5, 0.5}, {1, 1}, {0, 1} };

    for (auto p : index_range(points))
      mesh.add_point(points[p], /*id=*/p);

    std::unique_ptr<Elem> polygon =
        std::make_unique<C0Polygon>(cast_int<unsigned int>(points.size()));
    for (auto i : index_range(points))
      polygon->set_node(i, mesh.node_ptr(i));

    polygon->set_id() = 0;
    Elem *elem = mesh.add_elem(std::move(polygon));
    elem->subdomain_id() = 1;
    mesh.prepare_for_use();
  }

  void test_write_and_read_pentagon()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    this->build_pentagon(mesh);

    {
      ExodusII_IO exii(mesh);
      exii.write("write_exodus_C0POLYGON.e");
    }

    Mesh input_mesh(*TestCommWorld);
    ExodusII_IO exii_input(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii_input.read("write_exodus_C0POLYGON.e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(cast_int<dof_id_type>(1), input_mesh.n_elem());

    const Elem *elem = input_mesh.elem_ptr(0);
    CPPUNIT_ASSERT(elem);
    CPPUNIT_ASSERT_EQUAL(C0POLYGON, elem->type());
    CPPUNIT_ASSERT_EQUAL(5u, elem->n_nodes());

    for (auto i : make_range(5))
      CPPUNIT_ASSERT_EQUAL(cast_int<dof_id_type>(i), elem->node_id(i));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ExodusC0PolygonTest);
#endif // LIBMESH_DIM > 1

#if LIBMESH_DIM > 2
class ExodusC0PolyhedronTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(ExodusC0PolyhedronTest);

  CPPUNIT_TEST(test_write_cube_header);
  CPPUNIT_TEST(test_write_hexagonal_prism_header);
  CPPUNIT_TEST(test_write_and_read_hexagonal_prism);

  CPPUNIT_TEST_SUITE_END();

  ExodusHeaderInfo write_and_read_header(Mesh &mesh,
                                         const std::string &filename)
  {
    {
      ExodusII_IO exii(mesh);
      exii.write(filename);
    }

    TestCommWorld->barrier();

    Mesh header_mesh(*TestCommWorld);
    ExodusII_IO exii(header_mesh);
    return exii.read_header(filename);
  }

  void build_c0polyhedron(Mesh &mesh,
                          const std::vector<Point> &points,
                          const std::vector<std::vector<unsigned int>> &nodes_on_side)
  {
    for (auto p : index_range(points))
      mesh.add_point(points[p], /*id=*/p);

    std::vector<std::shared_ptr<Polygon>> sides(nodes_on_side.size());
    for (auto s : index_range(nodes_on_side))
      {
        const auto &nodes_on_s = nodes_on_side[s];
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
    Elem *elem = mesh.add_elem(std::move(polyhedron));
    elem->subdomain_id() = 1;
    mesh.cache_elem_data();
    mesh.prepare_for_use();
  }

  void test_write_cube_header()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_cube(mesh, 2, 2, 2,
                                      -1., 1.,
                                      -1., 1.,
                                      -1., 1.,
                                      C0POLYHEDRON);
    mesh.get_boundary_info().clear();

    for (auto &elem : mesh.element_ptr_range())
      elem->subdomain_id() = 1;
    mesh.cache_elem_data();

    ExodusHeaderInfo header_info =
        this->write_and_read_header(mesh, "write_exodus_C0POLYHEDRON.e");

    CPPUNIT_ASSERT_EQUAL(header_info.num_dim, 3);
    CPPUNIT_ASSERT_EQUAL(header_info.num_elem, 8);
    CPPUNIT_ASSERT_EQUAL(header_info.num_elem_blk, 1);
    CPPUNIT_ASSERT_EQUAL(header_info.num_face, 48);
    CPPUNIT_ASSERT_EQUAL(header_info.num_face_blk, 1);
    CPPUNIT_ASSERT_EQUAL(header_info.num_node_sets, 0);
    CPPUNIT_ASSERT_EQUAL(header_info.num_side_sets, 0);
  }

  void test_write_hexagonal_prism_header()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
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

    this->build_c0polyhedron(mesh, points, nodes_on_side);

    ExodusHeaderInfo header_info =
        this->write_and_read_header(mesh, "write_exodus_C0POLYHEDRON_HEXPRISM.e");

    CPPUNIT_ASSERT_EQUAL(header_info.num_dim, 3);
    CPPUNIT_ASSERT_EQUAL(header_info.num_elem, 1);
    CPPUNIT_ASSERT_EQUAL(header_info.num_elem_blk, 1);
    CPPUNIT_ASSERT_EQUAL(header_info.num_face, 8);
    CPPUNIT_ASSERT_EQUAL(header_info.num_face_blk, 1);
    CPPUNIT_ASSERT_EQUAL(header_info.num_node_sets, 0);
    CPPUNIT_ASSERT_EQUAL(header_info.num_side_sets, 0);
  }

  void test_write_and_read_hexagonal_prism()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
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

    this->build_c0polyhedron(mesh, points, nodes_on_side);

    {
      ExodusII_IO exii(mesh);
      exii.write("write_exodus_C0POLYHEDRON_HEXPRISM_READ.e");
    }

    Mesh input_mesh(*TestCommWorld);
    ExodusII_IO exii_input(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii_input.read("write_exodus_C0POLYHEDRON_HEXPRISM_READ.e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(cast_int<dof_id_type>(1), input_mesh.n_elem());

    const Elem *elem = input_mesh.elem_ptr(0);
    CPPUNIT_ASSERT(elem);
    CPPUNIT_ASSERT_EQUAL(C0POLYHEDRON, elem->type());
    CPPUNIT_ASSERT_EQUAL(12u, elem->n_vertices());
    CPPUNIT_ASSERT_EQUAL(8u, elem->n_sides());

    for (auto s : index_range(nodes_on_side))
      {
        const auto side_nodes = elem->nodes_on_side(s);
        CPPUNIT_ASSERT_EQUAL(nodes_on_side[s].size(), side_nodes.size());
        for (auto n : index_range(nodes_on_side[s]))
          CPPUNIT_ASSERT_EQUAL(cast_int<dof_id_type>(nodes_on_side[s][n]),
                               elem->node_id(side_nodes[n]));
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ExodusC0PolyhedronTest);

INSTANTIATE_EXODUSTEST(TET4);
INSTANTIATE_EXODUSTEST(TET10);
INSTANTIATE_EXODUSTEST(TET14);

INSTANTIATE_EXODUSTEST(HEX8);
INSTANTIATE_EXODUSTEST(HEX20);
INSTANTIATE_EXODUSTEST(HEX27);

INSTANTIATE_EXODUSTEST(PRISM6);
INSTANTIATE_EXODUSTEST(PRISM15);
INSTANTIATE_EXODUSTEST(PRISM18);
INSTANTIATE_EXODUSTEST(PRISM20);
INSTANTIATE_EXODUSTEST(PRISM21);

// These tests use PointLocator, which uses contains_point(), which
// uses inverse_map(), which doesn't play nicely on Pyramids unless we
// have exceptions support
#ifdef LIBMESH_ENABLE_EXCEPTIONS
INSTANTIATE_EXODUSTEST(PYRAMID5);
INSTANTIATE_EXODUSTEST(PYRAMID13);
INSTANTIATE_EXODUSTEST(PYRAMID14);
INSTANTIATE_EXODUSTEST(PYRAMID18);
#endif
#endif // LIBMESH_DIM > 2

#endif // LIBMESH_HAVE_EXODUS_API
