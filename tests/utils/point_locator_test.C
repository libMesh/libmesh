#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/node.h>
#include <libmesh/parallel.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;



class PointLocatorTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( PointLocatorTest );

  CPPUNIT_TEST( testLocatorOnEdge3 );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testLocatorOnQuad9 );
  CPPUNIT_TEST( testLocatorOnTri6 );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testLocatorOnHex27 );
  CPPUNIT_TEST( testPlanar );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testLocator(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    const unsigned int n_elem_per_side = 5;
    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    const unsigned int ymax = test_elem->dim() > 1;
    const unsigned int zmax = test_elem->dim() > 2;
    const unsigned int ny = ymax * n_elem_per_side;
    const unsigned int nz = zmax * n_elem_per_side;

    MeshTools::Generation::build_cube (mesh,
                                       n_elem_per_side,
                                       ny,
                                       nz,
                                       0., 1.,
                                       0., ymax,
                                       0., zmax,
                                       elem_type);

    std::unique_ptr<PointLocatorBase> locator = mesh.sub_point_locator();

    if (!mesh.is_serial())
      locator->enable_out_of_mesh_mode();

    for (unsigned int i=0; i != n_elem_per_side+1; ++i)
      {
        for (unsigned int j=0; j != ny+1; ++j)
          {
            for (unsigned int k=0; k != nz+1; ++k)
              {
                const libMesh::Real h = libMesh::Real(1)/n_elem_per_side;
                Point p(i*h, j*h, k*h);

                const Elem *elem = locator->operator()(p);

                bool found_elem = elem;
                if (!mesh.is_serial())
                  mesh.comm().max(found_elem);

                CPPUNIT_ASSERT(found_elem);
                if (elem)
                  {
                    CPPUNIT_ASSERT(elem->contains_point(p));
                  }

                const Node *node = locator->locate_node(p);

                bool found_node = node;
                if (!mesh.is_serial())
                  mesh.comm().max(found_node);

                CPPUNIT_ASSERT(found_node);

                if (node)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(i*h, (*node)(0),
                                            TOLERANCE*TOLERANCE);
                    if (LIBMESH_DIM > 1)
                      LIBMESH_ASSERT_FP_EQUAL(j*h, (*node)(1),
                                              TOLERANCE*TOLERANCE);
                    if (LIBMESH_DIM > 2)
                      LIBMESH_ASSERT_FP_EQUAL(k*h, (*node)(2),
                                              TOLERANCE*TOLERANCE);
                  }
              }
          }
      }
  }

  void testPlanar()
  {
    // Here we test locating points in a Mesh which lies slightly above the z-axis
    Mesh mesh(*TestCommWorld);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/10, /*ny=*/10,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        TRI3);

    // Move all nodes a small amount in the +z-direction
    for (auto & node : mesh.node_ptr_range())
      (*node)(2) += 3.1e-15;

    // Construct a PointLocator object
    std::unique_ptr<PointLocatorBase> locator = mesh.sub_point_locator();

    // Turn on out-of-mesh-mode to handle parallel testing
    if (!mesh.is_serial())
      locator->enable_out_of_mesh_mode();

    // Test locating a Point which is in the z=0 plane, i.e. slightly
    // below the plane of the Mesh, but which should otherwise be
    // found within the Mesh.
    Point p(0.53, 0.7, 0.0);
    const Elem *elem = (*locator)(p);

    bool found_elem = elem;
    if (!mesh.is_serial())
      mesh.comm().max(found_elem);

    CPPUNIT_ASSERT(found_elem);

    // Note: Tri3::contains_point() returns true for points which are
    // out-of-plane by less than TOLERANCE by default.
    if (elem)
      CPPUNIT_ASSERT(elem->contains_point(p));
  }

  void testLocatorOnEdge3() { testLocator(EDGE3); }
  void testLocatorOnQuad9() { testLocator(QUAD9); }
  void testLocatorOnTri6()  { testLocator(TRI6); }
  void testLocatorOnHex27() { testLocator(HEX27); }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PointLocatorTest );
