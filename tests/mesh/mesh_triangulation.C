#include <libmesh/parallel_implementation.h> // max()
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_triangle_holes.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/point.h>
#include <libmesh/poly2tri_triangulator.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>


using namespace libMesh;

class MeshTriangulationTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * interfaces to triangulation libraries
   */
public:
  CPPUNIT_TEST_SUITE( MeshTriangulationTest );

  CPPUNIT_TEST( testTriangleHoleArea );

#ifdef LIBMESH_HAVE_POLY2TRI
  CPPUNIT_TEST( testPoly2Tri );
  CPPUNIT_TEST( testPoly2TriHoles );
  CPPUNIT_TEST( testPoly2TriSegments );
  CPPUNIT_TEST( testPoly2TriRefined );
  CPPUNIT_TEST( testPoly2TriHolesRefined );
#endif

#ifdef LIBMESH_HAVE_TRIANGLE
  CPPUNIT_TEST( testTriangle );
  CPPUNIT_TEST( testTriangleHoles );
  CPPUNIT_TEST( testTriangleSegments );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testTriangleHoleArea()
  {
    // Using center=(1,0), radius=2 for the heck of it
    Point center{1};
    Real radius = 2;
    std::vector<TriangulatorInterface::PolygonHole> polyholes;

    // Line
    polyholes.emplace_back(center, radius, 2);
    // Triangle
    polyholes.emplace_back(center, radius, 3);
    // Square
    polyholes.emplace_back(center, radius, 4);
    // Pentagon
    polyholes.emplace_back(center, radius, 5);
    // Hexagon
    polyholes.emplace_back(center, radius, 6);

    for (int i=0; i != 5; ++i)
      {
        const int n_sides = i+2;
        const TriangulatorInterface::Hole & hole = polyholes[i];

        // Really?  This isn't until C++20?
        constexpr double my_pi = 3.141592653589793238462643383279;

        const Real computed_area = hole.area();
        const Real theta = my_pi/n_sides;
        const Real half_side_length = radius*std::cos(theta);
        const Real apothem = radius*std::sin(theta);
        const Real area = n_sides * apothem * half_side_length;

        LIBMESH_ASSERT_FP_EQUAL(computed_area, area, TOLERANCE*TOLERANCE);
      }

    TriangulatorInterface::ArbitraryHole arbhole {center, {{0,-1},{2,-1},{2,1},{0,2}}};
    LIBMESH_ASSERT_FP_EQUAL(arbhole.area(), Real(5), TOLERANCE*TOLERANCE);

#ifdef LIBMESH_HAVE_TRIANGLE
    // Make sure we're compatible with the old naming structure too
    TriangleInterface::PolygonHole square(center, radius, 4);
    LIBMESH_ASSERT_FP_EQUAL(square.area(), 2*radius*radius, TOLERANCE*TOLERANCE);
#endif
  }

  void testTriangulatorBase(MeshBase & mesh,
                            TriangulatorInterface & triangulator)
  {
    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    // Don't try to insert points yet
    triangulator.desired_area() = 1000;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(2));
    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        // Make sure we're not getting any inverted elements
        auto cross_prod =
          (elem->point(1) - elem->point(0)).cross
            (elem->point(2) - elem->point(0));

        CPPUNIT_ASSERT_GREATER(Real(0), cross_prod(2));

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
  }


  void testTriangulator(MeshBase & mesh,
                        TriangulatorInterface & triangulator)
  {
    // A non-square quad, so we don't have ambiguity about which
    // diagonal a Delaunay algorithm will pick.
    // Manually-numbered points, so we can use the point numbering as
    // a segment ordering even on DistributedMesh.
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,0), 1);
    mesh.add_point(Point(1,2), 2);
    mesh.add_point(Point(0,1), 3);

    this->testTriangulatorBase(mesh, triangulator);
  }



  void testTriangulatorHoles(MeshBase & mesh,
                             TriangulatorInterface & triangulator)
  {
    // A square quad; we'll put a diamond hole in the middle to make
    // the Delaunay selection unambiguous.
    mesh.add_point(Point(-1,-1), 0);
    mesh.add_point(Point(1,-1), 1);
    mesh.add_point(Point(1,1), 2);
    mesh.add_point(Point(-1,1), 3);

    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    // Add a diamond hole in the center
    TriangulatorInterface::PolygonHole diamond(Point(0), std::sqrt(2)/2, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };
    triangulator.attach_hole_list(&holes);

    // Don't try to insert points yet
    triangulator.desired_area() = 1000;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(8));

    // Center coordinates for all the elements we expect
    Real r2p2o6 = (std::sqrt(Real(2))+2)/6;
    Real r2p4o6 = (std::sqrt(Real(2))+4)/6;

    std::vector <Point> expected_centers
    { {r2p2o6,r2p2o6}, {r2p2o6,-r2p2o6},
      {-r2p2o6,r2p2o6}, {-r2p2o6,-r2p2o6},
      {0,r2p4o6}, {r2p4o6, 0},
      {0,-r2p4o6}, {-r2p4o6, 0}
    };

    std::vector<bool> found_centers(expected_centers.size(), false);

    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        // Make sure we're not getting any inverted elements
        auto cross_prod =
          (elem->point(1) - elem->point(0)).cross
            (elem->point(2) - elem->point(0));

        CPPUNIT_ASSERT_GREATER(Real(0), cross_prod(2));

        // Make sure we're finding all the elements we expect
        Point center = elem->vertex_average();

        bool found_mine = false;
        for (auto i : index_range(expected_centers))
          {
            Point possible = expected_centers[i];

            if (possible.absolute_fuzzy_equals(center, TOLERANCE*TOLERANCE))
              {
                found_mine = true;
                found_centers[i] = true;
              }
          }
        CPPUNIT_ASSERT(found_mine);
      }

    mesh.comm().max(found_centers);

    for (auto found_it : found_centers)
      CPPUNIT_ASSERT(found_it);
  }


  void testTriangulatorSegments(MeshBase & mesh,
                                TriangulatorInterface & triangulator)
  {
    // The same quad as testTriangulator, but out of order
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,2), 1);
    mesh.add_point(Point(1,0), 2);
    mesh.add_point(Point(0,1), 3);

    // Segments to put them in order
    triangulator.segments = {{0,2},{2,1},{1,3},{3,0}};

    this->testTriangulatorBase(mesh, triangulator);
  }


#ifdef LIBMESH_HAVE_TRIANGLE
  void testTriangle()
  {
    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulator(mesh, triangle);
  }


  void testTriangleHoles()
  {
    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorHoles(mesh, triangle);
  }


  void testTriangleSegments()
  {
    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorSegments(mesh, triangle);
  }
#endif // LIBMESH_HAVE_TRIANGLE


#ifdef LIBMESH_HAVE_POLY2TRI
  void testPoly2Tri()
  {
    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulator(mesh, p2t_tri);
  }


  void testPoly2TriHoles()
  {
    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorHoles(mesh, p2t_tri);
  }


  void testPoly2TriSegments()
  {
    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorSegments(mesh, p2t_tri);
  }


  void testPoly2TriRefinementBase
    (const std::vector<TriangulatorInterface::Hole*> * holes,
     Real expected_total_area,
     dof_id_type n_original_elem,
     Real desired_area = 0.1)
  {
    Mesh mesh(*TestCommWorld);
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,0), 1);
    mesh.add_point(Point(1,2), 2);
    mesh.add_point(Point(0,1), 3);

    Poly2TriTriangulator triangulator(mesh);

    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    if (holes)
      triangulator.attach_hole_list(holes);

    // Try to insert points!
    triangulator.desired_area() = desired_area;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

    CPPUNIT_ASSERT(mesh.n_elem() > n_original_elem);

    Real area = 0;
    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->level(), 0u);
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        const Real my_area = elem->volume();

        // my_area <= desired_area, wow this macro ordering hurts
        CPPUNIT_ASSERT_LESSEQUAL(desired_area, my_area);

        area += my_area;
      }

    mesh.comm().sum(area);

    LIBMESH_ASSERT_FP_EQUAL(area, expected_total_area, TOLERANCE*TOLERANCE);
  }

  void testPoly2TriRefined()
  {
    testPoly2TriRefinementBase(nullptr, 1.5, 2);
  }

  void testPoly2TriHolesRefined()
  {
    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    testPoly2TriRefinementBase(&holes, 1.25, 8);
  }
#endif // LIBMESH_HAVE_POLY2TRI

};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshTriangulationTest );
