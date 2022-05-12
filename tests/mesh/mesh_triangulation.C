#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_holes.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/parallel_implementation.h> // max()
#include <libmesh/parsed_function.h>
#include <libmesh/point.h>
#include <libmesh/poly2tri_triangulator.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <algorithm>
#include <cmath>


using namespace libMesh;

class MeshTriangulationTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * interfaces to triangulation libraries
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshTriangulationTest );

  CPPUNIT_TEST( testTriangleHoleArea );

#ifdef LIBMESH_HAVE_POLY2TRI
  CPPUNIT_TEST( testPoly2Tri );
  CPPUNIT_TEST( testPoly2TriInterp );
  CPPUNIT_TEST( testPoly2TriInterp3 );
  CPPUNIT_TEST( testPoly2TriHoles );
  CPPUNIT_TEST( testPoly2TriMeshedHoles );
  CPPUNIT_TEST( testPoly2TriEdges );
  CPPUNIT_TEST( testPoly2TriEdgesRefined );
  CPPUNIT_TEST( testPoly2TriSegments );
  CPPUNIT_TEST( testPoly2TriRefined );
  CPPUNIT_TEST( testPoly2TriExtraRefined );
  CPPUNIT_TEST( testPoly2TriHolesRefined );
  CPPUNIT_TEST( testPoly2TriHolesInteriorRefined );
  CPPUNIT_TEST( testPoly2TriHolesInteriorExtraRefined );
  // This covers an old poly2tri collinearity-tolerance bug
  CPPUNIT_TEST( testPoly2TriHolesExtraRefined );

  CPPUNIT_TEST( testPoly2TriNonUniformRefined );
  CPPUNIT_TEST( testPoly2TriHolesNonUniformRefined );
#endif

#ifdef LIBMESH_HAVE_TRIANGLE
  CPPUNIT_TEST( testTriangle );
  CPPUNIT_TEST( testTriangleInterp );
  CPPUNIT_TEST( testTriangleInterp3 );
  CPPUNIT_TEST( testTriangleHoles );
  CPPUNIT_TEST( testTriangleMeshedHoles );
  CPPUNIT_TEST( testTriangleSegments );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testTriangleHoleArea()
  {
    LOG_UNIT_TEST;

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



  void testTriangulatorInterp(MeshBase & mesh,
                              TriangulatorInterface & triangulator,
                              int interpolate_boundary_points,
                              dof_id_type n_expected_elem)
  {
    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    // A non-square quad, so we don't have ambiguity about which
    // diagonal a Delaunay algorithm will pick.
    // Manually-numbered points, so we can use the point numbering as
    // a segment ordering even on DistributedMesh.
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,0), 1);
    mesh.add_point(Point(1,2), 2);
    mesh.add_point(Point(0,1), 3);

    const Real expected_total_area = 1.5;

    // Interpolate points!
    triangulator.set_interpolate_boundary_points(interpolate_boundary_points);

    // Try to insert points!
    triangulator.desired_area() = 1000;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_expected_elem);

    Real area = 0;
    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->level(), 0u);
        CPPUNIT_ASSERT_EQUAL(elem->type(), TRI3);

        area += elem->volume();
      }

    mesh.comm().sum(area);

    LIBMESH_ASSERT_FP_EQUAL(area, expected_total_area, TOLERANCE*TOLERANCE);
  }


  void testFoundCenters(const MeshBase & mesh,
                        const std::vector<Point> & expected_centers)
  {
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

    testFoundCenters(mesh, expected_centers);
  }


  void testTriangulatorMeshedHoles(MeshBase & mesh,
                                   TriangulatorInterface & triangulator)
  {
    // A square quad; we'll put a square hole in the middle.  Offset
    // this to catch a potential bug I missed on the first try.
    mesh.add_point(Point(19,19), 0);
    mesh.add_point(Point(21,19), 1);
    mesh.add_point(Point(21,21), 2);
    mesh.add_point(Point(19,21), 3);

    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    // Add a square meshed hole in the center
    Mesh centermesh { mesh.comm() };
    MeshTools::Generation::build_square (centermesh, 2, 2, 19.5, 20.5, 19.5, 20.5, QUAD4);

    TriangulatorInterface::MeshedHole centerhole { centermesh };

    CPPUNIT_ASSERT_EQUAL(centerhole.n_points(), 8u);
    CPPUNIT_ASSERT_EQUAL(centerhole.area(), Real(1));
    Point inside = centerhole.inside();
    CPPUNIT_ASSERT_EQUAL(inside(0), Real(20));
    CPPUNIT_ASSERT_EQUAL(inside(1), Real(20));

    const std::vector<TriangulatorInterface::Hole*> holes { &centerhole };
    triangulator.attach_hole_list(&holes);

    // Don't try to insert points yet
    triangulator.desired_area() = 1000;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(12));

    // Center coordinates for all the elements we expect
    std::vector <Point> expected_centers
      { {-0.5, Real(2)/3}, {0, Real(5)/6},
        {0.5, Real(2)/3},
        {Real(2)/3, -0.5}, {Real(5)/6, 0},
        {Real(2)/3, 0.5},
        {-0.5, -Real(2)/3}, {0, -Real(5)/6},
        {0.5, -Real(2)/3},
        {-Real(2)/3, -0.5}, {-Real(5)/6, 0},
        {-Real(2)/3, 0.5} };

    // With the same offset
    for (auto & p : expected_centers)
      p += Point(20,20);

    testFoundCenters(mesh, expected_centers);
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
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulator(mesh, triangle);
  }


  void testTriangleInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorInterp(mesh, triangle, 2, 6);
  }


  void testTriangleInterp3()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorInterp(mesh, triangle, 3, 10);
  }


  void testTriangleHoles()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorHoles(mesh, triangle);
  }


  void testTriangleMeshedHoles()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorMeshedHoles(mesh, triangle);
  }


  void testTriangleSegments()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorSegments(mesh, triangle);
  }
#endif // LIBMESH_HAVE_TRIANGLE


#ifdef LIBMESH_HAVE_POLY2TRI
  void testPoly2Tri()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulator(mesh, p2t_tri);
  }


  void testPoly2TriInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorInterp(mesh, p2t_tri, 2, 6);
  }


  void testPoly2TriInterp3()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorInterp(mesh, p2t_tri, 3, 10);
  }


  void testPoly2TriHoles()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorHoles(mesh, p2t_tri);
  }


  void testPoly2TriMeshedHoles()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorMeshedHoles(mesh, p2t_tri);
  }


  void testPoly2TriEdgesMesh(MeshBase & mesh)
  {
    // The same quad as testTriangulator, but out of order
    auto node0 = mesh.add_point(Point(0,0), 0);
    auto node1 = mesh.add_point(Point(1,2), 1);
    auto node2 = mesh.add_point(Point(1,0), 2);
    auto node3 = mesh.add_point(Point(0,1), 3);

    // Edges, also out of order, but enough to put them in order
    auto edge13 = mesh.add_elem(Elem::build(EDGE2));
    edge13->set_node(0) = node1;
    edge13->set_node(1) = node3;
    auto edge02 = mesh.add_elem(Elem::build(EDGE2));
    edge02->set_node(0) = node0;
    edge02->set_node(1) = node2;
    auto edge30 = mesh.add_elem(Elem::build(EDGE2));
    edge30->set_node(0) = node3;
    edge30->set_node(1) = node0;
    auto edge21 = mesh.add_elem(Elem::build(EDGE2));
    edge21->set_node(0) = node2;
    edge21->set_node(1) = node1;

    mesh.prepare_for_use();
  }


  void testPoly2TriEdges()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);
    testPoly2TriEdgesMesh(mesh);

    this->testTriangulatorBase(mesh, triangulator);
  }


  void testPoly2TriEdgesRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);
    testPoly2TriEdgesMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 14);
  }


  void testPoly2TriSegments()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorSegments(mesh, p2t_tri);
  }


  void testPoly2TriTrapMesh (UnstructuredMesh & mesh)
  {
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,0), 1);
    mesh.add_point(Point(1,2), 2);
    mesh.add_point(Point(0,1), 3);
  }

  void testPoly2TriRefinementBase
    (UnstructuredMesh & mesh,
     const std::vector<TriangulatorInterface::Hole*> * holes,
     Real expected_total_area,
     dof_id_type n_original_elem,
     Real desired_area = 0.1,
     FunctionBase<Real> * area_func = nullptr)
  {
    Poly2TriTriangulator triangulator(mesh);

    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    if (holes)
      triangulator.attach_hole_list(holes);

    // Try to insert points!
    triangulator.desired_area() = desired_area;
    triangulator.set_desired_area_function(area_func);
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

  void testPoly2TriRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 15);
  }

  void testPoly2TriExtraRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 150, 0.01);
  }

  void testPoly2TriNonUniformRefined()
  {
    ParsedFunction<Real> var_area {"0.002*(1+2*x)*(1+2*y)"};
    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 150, 0, &var_area);
  }

  void testPoly2TriHolesRefined()
  {
    LOG_UNIT_TEST;

    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 13);
  }

  void testPoly2TriHolesInteriorRefinedBase
    (dof_id_type n_original_elem,
     Real desired_area)
  {
    // Add a diamond hole, disallowing refinement of it
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);

    CPPUNIT_ASSERT_EQUAL(diamond.refine_boundary_allowed(), true);

    diamond.set_refine_boundary_allowed(false);

    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    // Doing extra refinement here to ensure that we had the
    // *opportunity* to refine the hole boundaries.
    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, n_original_elem, desired_area);

    // Checking that we have more outer boundary sides than we started
    // with, and exactly the 4 hole boundary sides we started with.
    auto side_bcs = mesh.get_boundary_info().build_side_list();

    int n_outer_sides = std::count_if(side_bcs.begin(), side_bcs.end(),
                                      [](auto t){return std::get<2>(t) == 0;});
    CPPUNIT_ASSERT_GREATER(4, n_outer_sides); // n_outer_sides > 4
    int n_hole_sides = std::count_if(side_bcs.begin(), side_bcs.end(),
                                     [](auto t){return std::get<2>(t) == 1;});
    CPPUNIT_ASSERT_EQUAL(n_hole_sides, 4);
  }


  void testPoly2TriHolesInteriorRefined()
  {
    LOG_UNIT_TEST;
    testPoly2TriHolesInteriorRefinedBase(25, 0.05);
  }


  void testPoly2TriHolesInteriorExtraRefined()
  {
    LOG_UNIT_TEST;
    // 0.01 creates slivers triggering a poly2tri exception for me - RHS
    testPoly2TriHolesInteriorRefinedBase(60, 0.02);
  }


  void testPoly2TriHolesExtraRefined()
  {
    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 125, 0.01);
  }

  void testPoly2TriHolesNonUniformRefined()
  {
    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    ParsedFunction<Real> var_area {"0.002*(0.25+2*x)*(0.25+2*y)"};
    Mesh mesh(*TestCommWorld);
    testPoly2TriTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 150, 0, &var_area);
  }


#endif // LIBMESH_HAVE_POLY2TRI

};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshTriangulationTest );
