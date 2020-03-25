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
  CPPUNIT_TEST_SUITE( MeshGenerationTest );

  CPPUNIT_TEST( buildLineEdge2 );
  CPPUNIT_TEST( buildLineEdge3 );
  CPPUNIT_TEST( buildLineEdge4 );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( buildSquareTri3 );
  CPPUNIT_TEST( buildSquareTri6 );
  CPPUNIT_TEST( buildSquareQuad4 );
  CPPUNIT_TEST( buildSquareQuad8 );
  CPPUNIT_TEST( buildSquareQuad9 );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( buildCubeTet4 );
  CPPUNIT_TEST( buildCubeTet10 );
  CPPUNIT_TEST( buildCubeHex8 );
  CPPUNIT_TEST( buildCubeHex20 );
  CPPUNIT_TEST( buildCubeHex27 );
  CPPUNIT_TEST( buildCubePrism6 );
  CPPUNIT_TEST( buildCubePrism15 );
  CPPUNIT_TEST( buildCubePrism18 );
  CPPUNIT_TEST( buildCubePyramid5 );
  CPPUNIT_TEST( buildCubePyramid13 );
  CPPUNIT_TEST( buildCubePyramid14 );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testBuildLine(UnstructuredMesh & mesh, unsigned int n, ElemType type)
  {
    MeshTools::Generation::build_line (mesh, n, -1.0, 2.0, type);
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n);
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                         cast_int<dof_id_type>((Elem::type_to_n_nodes_map[type]-1)*n + 1));

    BoundingBox bbox = MeshTools::create_bounding_box(mesh);
    CPPUNIT_ASSERT_EQUAL(bbox.min()(0), Real(-1.0));
    CPPUNIT_ASSERT_EQUAL(bbox.max()(0), Real(2.0));
  }

  void testBuildSquare(UnstructuredMesh & mesh, unsigned int n, ElemType type)
  {
    MeshTools::Generation::build_square (mesh, n, n, -2.0, 3.0, -4.0, 5.0, type);
    if (Elem::type_to_n_sides_map[type] == 4)
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n));
    else
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), cast_int<dof_id_type>(n*n*2));

    const unsigned int n_nodes = Elem::type_to_n_nodes_map[type];

    switch (n_nodes)
      {
      case 3: // First-order elements
      case 4:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((n+1)*(n+1)));
        break;
      case 6: // Second-order elements
      case 9:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)));
        break;
      default:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1) - n*n));
      }

    // Our bounding boxes can be loose on higher order elements, but
    // we can at least assert that they're not too tight
    BoundingBox bbox = MeshTools::create_bounding_box(mesh);
    CPPUNIT_ASSERT(bbox.min()(0) <= Real(-2.0));
    CPPUNIT_ASSERT(bbox.max()(0) >= Real(3.0));
    CPPUNIT_ASSERT(bbox.min()(1) <= Real(-4.0));
    CPPUNIT_ASSERT(bbox.max()(1) >= Real(5.0));
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
        if (type == PRISM6 || type == PRISM15 || type == PRISM18)
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
      case 18: // Second-order prisms and hexes use the same nodes
      case 27:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1)));
        break;
      case 20:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) - n*n*n - 3*(n+1)*n*n));
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
      case 14:
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                             cast_int<dof_id_type>((2*n+1)*(2*n+1)*(2*n+1) + 8*n*n*n));
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
  }

  void buildLineEdge2 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildLine(rep, 6, EDGE2);
    testBuildLine(dist, 6, EDGE2);
  }

  void buildLineEdge3 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildLine(rep, 6, EDGE3);
    testBuildLine(dist, 6, EDGE3);
  }

  void buildLineEdge4 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildLine(rep, 6, EDGE3);
    testBuildLine(dist, 6, EDGE3);
  }

  void buildSquareTri3 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildSquare(rep, 4, TRI3);
    testBuildSquare(dist, 4, TRI3);
  }

  void buildSquareTri6 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildSquare(rep, 4, TRI6);
    testBuildSquare(dist, 4, TRI6);
  }

  void buildSquareQuad4 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildSquare(rep, 4, QUAD4);
    testBuildSquare(dist, 4, QUAD4);
  }

  void buildSquareQuad8 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildSquare(rep, 4, QUAD8);
    testBuildSquare(dist, 4, QUAD8);
  }

  void buildSquareQuad9 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildSquare(rep, 4, QUAD9);
    testBuildSquare(dist, 4, QUAD9);
  }

  void buildCubeTet4 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, TET4);
    testBuildCube(dist, 2, TET4);
  }

  void buildCubeTet10 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, TET10);
    testBuildCube(dist, 2, TET10);
  }

  void buildCubeHex8 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, HEX8);
    testBuildCube(dist, 2, HEX8);
  }

  void buildCubeHex20 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, HEX20);
    testBuildCube(dist, 2, HEX20);
  }

  void buildCubeHex27 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, HEX27);
    testBuildCube(dist, 2, HEX27);
  }

  void buildCubePrism6 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PRISM6);
    testBuildCube(dist, 2, PRISM6);
  }

  void buildCubePrism15 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PRISM15);
    testBuildCube(dist, 2, PRISM15);
  }

  void buildCubePrism18 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PRISM18);
    testBuildCube(dist, 2, PRISM18);
  }
 
  void buildCubePyramid5 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PYRAMID5);
    testBuildCube(dist, 2, PYRAMID5);
  }
   
  void buildCubePyramid13 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PYRAMID13);
    testBuildCube(dist, 2, PYRAMID13);
  }
 
  void buildCubePyramid14 () {
    ReplicatedMesh rep(*TestCommWorld);
    DistributedMesh dist(*TestCommWorld);
    testBuildCube(rep, 2, PYRAMID14);
    testBuildCube(dist, 2, PYRAMID14);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshGenerationTest );
