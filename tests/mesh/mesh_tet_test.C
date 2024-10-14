#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_netgen_interface.h>
#include <libmesh/mesh_tetgen_interface.h>
#include <libmesh/mesh_tet_interface.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/parallel_implementation.h> // max()

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <algorithm>
#include <cmath>
#include <regex>


using namespace libMesh;

namespace {

Real build_octahedron (UnstructuredMesh & mesh, bool flip_tris,
                       Real xmin, Real xmax,
                       Real ymin, Real ymax,
                       Real zmin, Real zmax)
{
  MeshTools::Generation::surface_octahedron
    (mesh, xmin, xmax, ymin, ymax, zmin, zmax);

  // Octahedron volume
  return (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/6;
}

}


class MeshTetTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * interfaces to tetrahedralization libraries
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshTetTest );

#ifdef LIBMESH_HAVE_NETGEN
  // The most basic test to start with
  CPPUNIT_TEST( testNetGen );
  CPPUNIT_TEST( testNetGenError );
  CPPUNIT_TEST( testNetGenTets );
  CPPUNIT_TEST( testNetGenFlippedTris );
  CPPUNIT_TEST( testNetGenHole );

#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testNetGenSphereShell );
#endif

  // We'll get to more advanced features later
  /*
  CPPUNIT_TEST( testNetGenInterp );
  CPPUNIT_TEST( testNetGenInterp2 );
  CPPUNIT_TEST( testNetGenRefined );
  CPPUNIT_TEST( testNetGenNonRefined );
  CPPUNIT_TEST( testNetGenExtraRefined );
  */
#endif

  // Still need to work out the basics here - non-convex domains,
  // precise control of refinement, etc.
#ifdef LIBMESH_HAVE_TETGEN
  /*
  CPPUNIT_TEST( testTetGen );
  CPPUNIT_TEST( testTetGenError );
  CPPUNIT_TEST( testTetGenInterp );
  CPPUNIT_TEST( testTetGenInterp2 );
  */
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testExceptionBase(const char * re,
                         MeshBase & mesh,
                         MeshTetInterface & tetinterface,
                         dof_id_type expected_n_elem = DofObject::invalid_id,
                         dof_id_type expected_n_nodes = DofObject::invalid_id,
                         Real expected_volume = 0)
  {
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // We can't just CPPUNIT_ASSERT_THROW, because we want to make
    // sure we were thrown from the right place with the right error
    // message!
    bool threw_desired_exception = false;
    try {
      this->testTetInterfaceBase(mesh, tetinterface, expected_n_elem,
                                 expected_n_nodes, expected_volume);
    }
    catch (libMesh::LogicError & e) {
      std::regex msg_regex(re);
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_desired_exception = true;
    }
    catch (CppUnit::Exception & e) {
      throw e;
    }
    catch (...) {
      CPPUNIT_ASSERT_MESSAGE("Unexpected exception type thrown", false);
    }
    CPPUNIT_ASSERT(threw_desired_exception);
#endif
  }


  void testTetInterfaceBase(MeshBase & mesh,
                            MeshTetInterface & triangulator,
                            dof_id_type expected_n_elem = DofObject::invalid_id,
                            dof_id_type expected_n_nodes = DofObject::invalid_id,
                            Real expected_volume = 0)
  {
    triangulator.triangulate();

    if (expected_n_elem != DofObject::invalid_id)
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), expected_n_elem);

    if (expected_n_nodes != DofObject::invalid_id)
      CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), expected_n_nodes);

    if (expected_volume != 0)
      LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(mesh),
                              expected_volume,
                              TOLERANCE*TOLERANCE);

    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), TET4);

        // Make sure we're not getting any inverted elements
        CPPUNIT_ASSERT(!elem->is_flipped());
      }
  }


  void testHole(UnstructuredMesh & mesh,
                MeshTetInterface & triangulator)
  {
    std::unique_ptr<UnstructuredMesh> holemesh =
      std::make_unique<Mesh>(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh, 1, 1, 1,
                                       -2, 2, -2, 2, -2, 2);

    const Real hole_volume =
      build_octahedron(*holemesh, false, -1, 1, -1, 1, -1, 1);

    auto holes =
      std::make_unique<std::vector<std::unique_ptr<UnstructuredMesh>>>();

    holes->push_back(std::move(holemesh));

    triangulator.attach_hole_list(std::move(holes));

    const Real expected_volume =
      MeshTools::volume(mesh) - hole_volume;
    this->testTetInterfaceBase(mesh, triangulator, 32, 14,
                               expected_volume);
  }


#ifdef LIBMESH_ENABLE_AMR
  void testSphereShell(UnstructuredMesh & mesh,
                       MeshTetInterface & triangulator)
  {
    std::unique_ptr<UnstructuredMesh> holemesh =
      std::make_unique<Mesh>(*TestCommWorld);

    MeshTools::Generation::build_sphere (*holemesh, 1, 2, TET4);

    MeshTools::Generation::build_sphere (mesh, 1.25, 2, TET4);

    auto holes =
      std::make_unique<std::vector<std::unique_ptr<UnstructuredMesh>>>();

    holes->push_back(std::move(holemesh));

    triangulator.attach_hole_list(std::move(holes));

    // Netgen can't seem to triangulate this without inserting points,
    // so let MeshNetgenInterface know that we're allowed to insert
    // points
    triangulator.desired_volume() = 1000;

    this->testTetInterfaceBase(mesh, triangulator);
  }
#endif


  void testTrisToTets(UnstructuredMesh & mesh,
                      MeshTetInterface & triangulator,
                      bool flip_tris = false)
  {
    // An asymmetric octahedron, so we hopefully have an unambiguous
    // choice of shortest diagonal for a Delaunay algorithm to pick.
    const Real expected_volume =
      build_octahedron(mesh, flip_tris, -1, 1, -1, 1, -0.1, 0.1);

    this->testTetInterfaceBase(mesh, triangulator, /* n_elem = */ 4,
                               /* n_nodes = */ 6, expected_volume);
  }


  void testTrisToTetsError(UnstructuredMesh & mesh,
                           MeshTetInterface & triangulator,
                           bool flip_tris = false)
  {
    const Real expected_volume =
      build_octahedron(mesh, flip_tris, -1, 1, -1, 1, -0.1, 0.1);

    // Remove one tri, breaking the mesh
    for (auto elem : mesh.element_ptr_range())
      {
        Point center = elem->vertex_average();
        if (center(0) > 0 &&
            center(1) > 0 &&
            center(2) > 0)
          mesh.delete_elem(elem);
      }
    mesh.prepare_for_use();

    this->testExceptionBase("element with a null neighbor", mesh, triangulator,
                            /* n_elem = */ 4, /* n_nodes = */ 6,
                            expected_volume);
  }


  void testTetsToTets(MeshBase & mesh,
                      MeshTetInterface & triangulator)
  {
    // An asymmetric octahedron, so we hopefully have an unambiguous
    // choice of shortest diagonal for a Delaunay algorithm to pick.
    mesh.add_point(Point(0,0,-0.1), 0);
    mesh.add_point(Point(1,0,0), 1);
    mesh.add_point(Point(0,1,0), 2);
    mesh.add_point(Point(-1,0,0), 3);
    mesh.add_point(Point(0,-1,0), 4);
    mesh.add_point(Point(0,0,0.1), 5);

    auto add_tet = [&mesh](std::array<dof_id_type,4> nodes)
    {
      auto elem = mesh.add_elem(Elem::build(TET4));
      elem->set_node(0) = mesh.node_ptr(nodes[0]);
      elem->set_node(1) = mesh.node_ptr(nodes[1]);
      elem->set_node(2) = mesh.node_ptr(nodes[2]);
      elem->set_node(3) = mesh.node_ptr(nodes[3]);
    };

    // Split along a different diagonal to start
    add_tet({1,3,4,5});
    add_tet({1,3,5,2});
    add_tet({1,3,2,0});
    add_tet({1,3,0,4});

    mesh.prepare_for_use();

    const Real expected_volume = MeshTools::volume(mesh);

    this->testTetInterfaceBase(mesh, triangulator, /* n_elem = */ 4,
                               /* n_nodes = */ 6, expected_volume);
  }


#ifdef LIBMESH_HAVE_TETGEN
  void testTetGen()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TetGenMeshInterface tet_tet(mesh);
    testTrisToTets(mesh, tet_tet);
  }


  void testTetGenError()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TetGenMeshInterface tet_tet(mesh);
      testTrisToTetsError(mesh, tet_tet);
  }



  /*
  void testTetGenInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TetGenMeshInterface tet_tet(mesh);
    testTrisToTetsInterp(mesh, tet_tet, 1, 6);
  }


  void testTetGenInterp2()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TetGenMeshInterface tet_tet(mesh);
    testTrisToTetsInterp(mesh, tet_tet, 2, 10);
  }
  */

#endif // LIBMESH_HAVE_TETGEN


#ifdef LIBMESH_HAVE_NETGEN
  void testNetGen()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTrisToTets(mesh, net_tet);
  }


  void testNetGenError()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTrisToTetsError(mesh, net_tet);
  }


 void testNetGenTets()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTetsToTets(mesh, net_tet);
  }


  void testNetGenFlippedTris()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTrisToTets(mesh, net_tet, true);
  }


  void testNetGenHole()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testHole(mesh, net_tet);
  }


#ifdef LIBMESH_ENABLE_AMR
  void testNetGenSphereShell()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testSphereShell(mesh, net_tet);
  }
#endif


  /*
  void testNetGenInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTrisToTetsInterp(mesh, net_tet, 1, 6);
  }


  void testNetGenInterp2()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    NetGenMeshInterface net_tet(mesh);
    testTrisToTetsInterp(mesh, net_tet, 2, 10);
  }


  void testNetGenRefinementBase
    (UnstructuredMesh & mesh,
     const std::vector<MeshTetInterface::Hole*> * holes,
     Real expected_total_area,
     dof_id_type n_original_elem,
     Real desired_area = 0.1,
     FunctionBase<Real> * area_func = nullptr)
  {
    NetGenMeshInterface triangulator(mesh);

    if (holes)
      triangulator.attach_hole_list(holes);

    // Try to insert points!
    triangulator.desired_area() = desired_area;
    triangulator.set_desired_area_function(area_func);

    triangulator.triangulate();

    // If refinement should have increased our element count, check it
    if (desired_area || area_func)
      CPPUNIT_ASSERT_GREATER(n_original_elem, mesh.n_elem()); // n_elem+++
    else
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_original_elem);

    Real area = 0;
    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->level(), 0u);
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        const Real my_area = elem->volume();

        // my_area <= desired_area, wow this macro ordering hurts
        if (desired_area != 0)
          CPPUNIT_ASSERT_LESSEQUAL(desired_area, my_area);

        if (area_func != nullptr)
          for (auto v : make_range(elem->n_vertices()))
            {
              const Real local_desired_area =
                (*area_func)(elem->point(v));
              CPPUNIT_ASSERT_LESSEQUAL(local_desired_area, my_area);
            }

        area += my_area;
      }

    mesh.comm().sum(area);

    LIBMESH_ASSERT_FP_EQUAL(area, expected_total_area, TOLERANCE*TOLERANCE);
  }

  void testNetGenRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 15);
  }

  void testNetGenNonRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    // Make sure we see 0 as "don't refine", not "infinitely refine"
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 2, 0);
  }


  void testNetGenExtraRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 150, 0.01);
  }
*/

#endif // LIBMESH_HAVE_NETGEN

};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshTetTest );
