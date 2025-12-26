#include <libmesh/libmesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshGenerationTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of
   * MeshGeneration functions, as well as to indirectly verify the
   * MeshBase functions they rely on.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshGenerationTest );

  CPPUNIT_TEST( buildLineEdge2 );
  CPPUNIT_TEST( buildLineEdge3 );
  CPPUNIT_TEST( buildLineEdge4 );
#  ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( buildSphereEdge2 );
  CPPUNIT_TEST( buildSphereEdge3 );
// CPPUNIT_TEST( buildSphereEdge4 ); Doesn't work with AMR yet
#  endif

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( buildSquareTri3 );
  CPPUNIT_TEST( buildSquareTri6 );
  CPPUNIT_TEST( buildSquareTri7 );
  CPPUNIT_TEST( buildSquareQuad4 );
  CPPUNIT_TEST( buildSquareQuad8 );
  CPPUNIT_TEST( buildSquareQuad9 );
  CPPUNIT_TEST( buildSquareC0PolygonEven );
  CPPUNIT_TEST( buildSquareC0PolygonOdd );
#  ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( buildSphereTri3 );
  CPPUNIT_TEST( buildSphereQuad4 );
#  endif
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( buildCubeTet4 );
  CPPUNIT_TEST( buildCubeTet10 );
  CPPUNIT_TEST( buildCubeTet14 );
  CPPUNIT_TEST( buildCubeHex8 );
  CPPUNIT_TEST( buildCubeHex20 );
  CPPUNIT_TEST( buildCubeHex27 );
  CPPUNIT_TEST( buildCubePrism6 );
  CPPUNIT_TEST( buildCubePrism15 );
  CPPUNIT_TEST( buildCubePrism18 );
  CPPUNIT_TEST( buildCubePrism20 );
  CPPUNIT_TEST( buildCubePrism21 );

  // These tests throw an exception from contains_point() calls, and
  // this simply aborts() when exceptions are not enabled.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  CPPUNIT_TEST( buildCubePyramid5 );
  CPPUNIT_TEST( buildCubePyramid13 );
  CPPUNIT_TEST( buildCubePyramid14 );

#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( buildSphereHex27 );
#endif
#endif

#  ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( buildSphereHex8 );
#  endif
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  std::unique_ptr<UnstructuredMesh> new_mesh (bool is_replicated)
  {
    if (is_replicated)
      return std::make_unique<ReplicatedMesh>(*TestCommWorld);
    return std::make_unique<DistributedMesh>(*TestCommWorld);
  }

public:
  void setUp() {}

  void tearDown() {}

  void testBuildLine(UnstructuredMesh & mesh, unsigned int n, ElemType type)
  {
    MeshTools::Generation::build_line (mesh, n, -1.0, 2.0, type);
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                         cast_int<dof_id_type>((Elem::type_to_n_nodes_map[type]-1)*n + 1));

    BoundingBox bbox = MeshTools::create_bounding_box(mesh);
    CPPUNIT_ASSERT_EQUAL(bbox.min()(0), Real(-1.0));
    CPPUNIT_ASSERT_EQUAL(bbox.max()(0), Real(2.0));

    // Do serial assertions *after* all parallel assertions, so we
    // stay in sync after failure on only some processor(s)
    for (auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT(elem->has_affine_map());
  }

  void testBuildSquare(UnstructuredMesh & mesh, unsigned int n, ElemType type)
  {
    MeshTools::Generation::build_square (mesh, n, n, -2.0, 3.0, -4.0, 5.0, type);
    if (type == C0POLYGON)
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n + 4 + 2 * (n - 1) + ((n - 1) / 2)));
    else if (Elem::type_to_n_sides_map[type] == 4)
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n));
    else
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*2));

    switch (type)
      {
      case TRI3: // First-order elements
      case QUAD4:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((n+1)*(n+1)));
        break;
      case TRI6: // Second-order elements
      case QUAD9:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)));
        break;
      case QUAD8: // Not-really-second-order element missing center nodes
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1) - n*n));
        break;
      case TRI7: // Not-really-second-order element with extra center nodes
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1) + 2*n*n));
        break;
      case C0POLYGON:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>(4 + 2*n*n + (n - 1) + 2*n + 2 * (n%2)));
        break;
      default: // Wait, what did we try to build?
        CPPUNIT_ASSERT(false);
      }

    // Our bounding boxes can be loose on higher order elements, but
    // we can at least assert that they're not too tight
    BoundingBox bbox = MeshTools::create_bounding_box(mesh);
    CPPUNIT_ASSERT(bbox.min()(0) <= Real(-2.0));
    CPPUNIT_ASSERT(bbox.max()(0) >= Real(3.0));
    CPPUNIT_ASSERT(bbox.min()(1) <= Real(-4.0));
    CPPUNIT_ASSERT(bbox.max()(1) >= Real(5.0));

    // Do serial assertions *after* all parallel assertions, so we
    // stay in sync after failure on only some processor(s)
    if (type != C0POLYGON)
      for (auto & elem : mesh.element_ptr_range())
        CPPUNIT_ASSERT(elem->has_affine_map());
  }

  void testBuildCube(UnstructuredMesh & mesh, unsigned int n, ElemType type)
  {
    MeshTools::Generation::build_cube (mesh, n, n, n, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, type);
    switch (Elem::type_to_n_sides_map[type])
      {
      case 4: // tets
        CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*n*24));
        break;
      case 5: // prisms, pyramids
        if (type == PRISM6 || type == PRISM15 || type == PRISM18 ||
            type == PRISM20 || type == PRISM21)
          CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*n*2));
        else
          CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*n*6));
        break;
      case 6: // hexes
        CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*n));
        break;
      default:
        libmesh_error();
      }


    switch (Elem::type_to_n_nodes_map[type])
      {
      case 4: // First-order tets
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((n+1)*(n+1)*(n+1) + n*n*n + 3*(n+1)*n*n));
        break;
      case 6: // First-order prisms and hexes use the same nodes
      case 8:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((n+1)*(n+1)*(n+1)));
        break;
      case 10: // Second-order tets
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 14*n*n*n + 4*3*(n+1)*n*n));
        break;
      case 18:
      case 27: // Second-order prisms and hexes use the same nodes
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1)));
        break;
      case 20:
        if (type == HEX20)
          CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                               cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) - n*n*n - 3*(n+1)*n*n));
        if (type == PRISM20)
          CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                               cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 2*(n+1)*n*n));
        break;
      case 21: // Prisms based on full Tri7 cross sections
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 2*(2*n+1)*n*n));
        break;
      case 15: // weird partial order prism
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) - n*n*n - 2*(n+1)*n*n));
        break;
      case 5: // pyramids
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((n+1)*(n+1)*(n+1) + n*n*n));
        break;
      case 13:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 8*n*n*n - 3*(n+1)*n*n));
        break;
      case 14: // pyramids, tets
        if (type == PYRAMID14)
          CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                               cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 8*n*n*n));
        else // TET14
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 14*n*n*n + 4*3*(n+1)*n*n +
                                                   36*n*n*n + 4*3*(n+1)*n*n));
        break;
      default:
        libmesh_error();
      }

    // Our bounding boxes can be loose on higher order elements, but
    // we can at least assert that they're not too tight
    BoundingBox bbox = MeshTools::create_bounding_box(mesh);
    CPPUNIT_ASSERT(bbox.min()(0) <= Real(-2.0));
    CPPUNIT_ASSERT(bbox.max()(0) >= Real(3.0));
    CPPUNIT_ASSERT(bbox.min()(1) <= Real(-4.0));
    CPPUNIT_ASSERT(bbox.max()(1) >= Real(5.0));
    CPPUNIT_ASSERT(bbox.min()(2) <= Real(-6.0));
    CPPUNIT_ASSERT(bbox.max()(2) >= Real(7.0));

    // We don't yet try to do affine map optimizations on pyramids
    if (type == PYRAMID5 ||
        type == PYRAMID13 ||
        type == PYRAMID14)
      return;

    // Do serial assertions *after* all parallel assertions, so we
    // stay in sync after failure on only some processor(s)
    for (auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT(elem->has_affine_map());
  }

  void testBuildSphere(unsigned int n_ref, ElemType type)
  {
    ReplicatedMesh rmesh(*TestCommWorld);
    MeshTools::Generation::build_sphere (rmesh, 2.0, n_ref, type);

    DistributedMesh dmesh(*TestCommWorld);
    dmesh.allow_renumbering(false);
    MeshTools::Generation::build_sphere (dmesh, 2.0, n_ref, type);
  }


  typedef void (MeshGenerationTest::*Builder)(UnstructuredMesh&, unsigned int, ElemType);

  void tester(Builder f, unsigned int n, ElemType type)
  {
    for (int is_replicated = 0; is_replicated != 2; ++is_replicated)
      {
        for (int skip_renumber = 0 ; skip_renumber != 2; ++skip_renumber)
          {
            std::unique_ptr<UnstructuredMesh> mesh =
              new_mesh(is_replicated);
            mesh->allow_renumbering(!skip_renumber);
            (this->*f)(*mesh, n, type);
          }
      }
  }

  void buildLineEdge2 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildLine, 5, EDGE2); }
  void buildLineEdge3 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildLine, 5, EDGE3); }
  void buildLineEdge4 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildLine, 5, EDGE4); }

  void buildSphereEdge2 ()     { LOG_UNIT_TEST; testBuildSphere(2, EDGE2); }
  void buildSphereEdge3 ()     { LOG_UNIT_TEST; testBuildSphere(2, EDGE3); }
  void buildSphereEdge4 ()     { LOG_UNIT_TEST; testBuildSphere(2, EDGE4); }

  void buildSquareTri3 ()    { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 3, TRI3); }
  void buildSquareTri6 ()    { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 4, TRI6); }
  void buildSquareTri7 ()    { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 4, TRI7); }
  void buildSquareQuad4 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 4, QUAD4); }
  void buildSquareQuad8 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 4, QUAD8); }
  void buildSquareQuad9 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 4, QUAD9); }
  void buildSquareC0PolygonOdd() { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 5, C0POLYGON); }
  void buildSquareC0PolygonEven(){ LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildSquare, 6, C0POLYGON); }

  void buildSphereTri3 ()     { LOG_UNIT_TEST; testBuildSphere(2, TRI3); }
  void buildSphereQuad4 ()     { LOG_UNIT_TEST; testBuildSphere(2, QUAD4); }

  void buildCubeTet4 ()      { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, TET4); }
  void buildCubeTet10 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, TET10); }
  void buildCubeTet14 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, TET14); }
  void buildCubeHex8 ()      { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, HEX8); }
  void buildCubeHex20 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, HEX20); }
  void buildCubeHex27 ()     { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, HEX27); }
  void buildCubePrism6 ()    { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PRISM6); }
  void buildCubePrism15 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PRISM15); }
  void buildCubePrism18 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PRISM18); }
  void buildCubePrism20 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PRISM20); }
  void buildCubePrism21 ()   { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PRISM21); }

  // These tests throw an exception from contains_point() calls, and
  // this simply aborts() when exceptions are not enabled.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  void buildCubePyramid5 ()  { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PYRAMID5); }
  void buildCubePyramid13 () { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PYRAMID13); }
  void buildCubePyramid14 () { LOG_UNIT_TEST; tester(&MeshGenerationTest::testBuildCube, 2, PYRAMID14); }
#endif

  void buildSphereHex8 ()     { LOG_UNIT_TEST; testBuildSphere(2, HEX8); }
  void buildSphereHex27 ()     { LOG_UNIT_TEST; testBuildSphere(2, HEX27); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshGenerationTest );
