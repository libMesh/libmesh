#include <libmesh/elem.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/mesh.h>
#include <libmesh/point.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshTriangulationTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * interfaces to triangulation libraries
   */
public:
  CPPUNIT_TEST_SUITE( MeshTriangulationTest );

#ifdef LIBMESH_HAVE_TRIANGLE
  CPPUNIT_TEST( testTriangle );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testTriangle()
  {
#ifdef LIBMESH_HAVE_TRIANGLE
    Mesh mesh(*TestCommWorld);

    // A non-square quad, so we don't have ambiguity about which
    // diagonal a Delaunay algorithm will pick.
    mesh.add_point(Point(0,0));
    mesh.add_point(Point(1,0));
    mesh.add_point(Point(1,2));
    mesh.add_point(Point(0,1));

    TriangleInterface triangle(mesh);

    // Don't try to insert points yet
    triangle.desired_area() = 1000;
    triangle.minimum_angle() = 0;
    triangle.smooth_after_generating() = false;

    triangle.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(2));
    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        bool found_triangle = false;
        for (const auto & node : elem->node_ref_range())
          {
            const Point & point = node;
            if (point == Point(0,0))
              {
                found_triangle = true;
                CPPUNIT_ASSERT((elem->point(0) == Point(0,0) &&
                                elem->point(1) == Point(1,0) &&
                                elem->point(2) == Point(0,1)) ||
                               (elem->point(1) == Point(0,0) &&
                                elem->point(2) == Point(1,0) &&
                                elem->point(0) == Point(0,1)) ||
                               (elem->point(2) == Point(0,0) &&
                                elem->point(0) == Point(1,0) &&
                                elem->point(1) == Point(0,1)));
              }
            if (point == Point(1,2))
              {
                found_triangle = true;
                CPPUNIT_ASSERT((elem->point(0) == Point(0,1) &&
                                elem->point(1) == Point(1,0) &&
                                elem->point(2) == Point(1,2)) ||
                               (elem->point(1) == Point(0,1) &&
                                elem->point(2) == Point(1,0) &&
                                elem->point(0) == Point(1,2)) ||
                               (elem->point(2) == Point(0,1) &&
                                elem->point(0) == Point(1,0) &&
                                elem->point(1) == Point(1,2)));
              }
          }
        CPPUNIT_ASSERT(found_triangle);
      }
#endif // LIBMESH_HAVE_TRIANGLE
  }

};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshTriangulationTest );
