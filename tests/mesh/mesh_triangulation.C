#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
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
#include <regex>


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
  CPPUNIT_TEST( testTriangleHoleContains );

#ifdef LIBMESH_HAVE_POLY2TRI
  CPPUNIT_TEST( testPoly2Tri );
  CPPUNIT_TEST( testPoly2TriHalfDomain );
  CPPUNIT_TEST( testPoly2TriHalfDomainEdge3 );
  CPPUNIT_TEST( testPoly2TriInterp );
  CPPUNIT_TEST( testPoly2TriInterp2 );
  CPPUNIT_TEST( testPoly2TriHoles );
  CPPUNIT_TEST( testPoly2TriMeshedHoles );
#  ifdef LIBMESH_ENABLE_AMR
    CPPUNIT_TEST( testPoly2TriRoundHole );
#  endif
  CPPUNIT_TEST( testPoly2TriEdges );
  CPPUNIT_TEST( testPoly2TriEdge3s );
  CPPUNIT_TEST( testPoly2TriBadEdges );
  CPPUNIT_TEST( testPoly2TriBad1DMultiBoundary );
  CPPUNIT_TEST( testPoly2TriBad2DMultiBoundary );
  CPPUNIT_TEST( testPoly2TriEdgesRefined );
  CPPUNIT_TEST( testPoly2TriSegments );
  CPPUNIT_TEST( testPoly2TriRefined );
  CPPUNIT_TEST( testPoly2TriNonRefined );
  CPPUNIT_TEST( testPoly2TriExtraRefined );
  CPPUNIT_TEST( testPoly2TriHolesRefined );
  CPPUNIT_TEST( testPoly2TriHolesInterpRefined );
  CPPUNIT_TEST( testPoly2TriHolesInteriorRefined );
  CPPUNIT_TEST( testPoly2TriHolesInteriorExtraRefined );
  // This covers an old poly2tri collinearity-tolerance bug
  CPPUNIT_TEST( testPoly2TriHolesExtraRefined );

  CPPUNIT_TEST( testPoly2TriNonUniformRefined );
  CPPUNIT_TEST( testPoly2TriHolesNonUniformRefined );
#endif

#ifdef LIBMESH_HAVE_TRIANGLE
  CPPUNIT_TEST( testTriangle );
  CPPUNIT_TEST( testTriangleHalfDomain );
  CPPUNIT_TEST( testTriangleInterp );
  CPPUNIT_TEST( testTriangleInterp2 );
  CPPUNIT_TEST( testTriangleHoles );
  CPPUNIT_TEST( testTriangleMeshedHoles );
#  ifdef LIBMESH_ENABLE_AMR
    CPPUNIT_TEST( testTriangleRoundHole );
#  endif
  CPPUNIT_TEST( testTriangleEdges );
  CPPUNIT_TEST( testTriangleSegments );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void commonSettings (TriangulatorInterface & triangulator)
  {
    // Use the point order to define the boundary, because our
    // Poly2Tri implementation doesn't do convex hulls yet, even when
    // that would give the same answer.
    triangulator.triangulation_type() = TriangulatorInterface::PSLG;

    // Don't try to insert points unless we're requested to later
    triangulator.desired_area() = 1000;
    triangulator.minimum_angle() = 0;
    triangulator.smooth_after_generating() = false;
    triangulator.set_verify_hole_boundaries(true);
  }

  void testExceptionBase(MeshBase & mesh,
                         TriangulatorInterface & triangulator,
                         const char * re)
  {
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // We can't just CPPUNIT_ASSERT_THROW, because we want to make
    // sure we were thrown from the right place with the right error
    // message!
    bool threw_desired_exception = false;
    try {
      this->testTriangulatorBase(mesh, triangulator);
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


  void testTriangleHoleContains()
  {
    LOG_UNIT_TEST;

    // Using center=(1,2), radius=2 for the heck of it
    Point center{1,2};
    Real radius = 2;
    std::vector<TriangulatorInterface::PolygonHole> polyholes;

    auto check_corners = [center, radius]
      (const TriangulatorInterface::Hole & hole,
       unsigned int np)
      {
        CPPUNIT_ASSERT_EQUAL(np, hole.n_points());
        for (auto i : make_range(np))
          {
            const Real theta = i * 2 * libMesh::pi / np;
            const Real xin = center(0) + radius * .99 * std::cos(theta);
            const Real xout = center(0) + radius * 1.01 * std::cos(theta);
            const Real yin = center(1) + radius * .99 * std::sin(theta);
            const Real yout = center(1) + radius * 1.01 * std::sin(theta);

            CPPUNIT_ASSERT(hole.contains(Point(xin, yin)));
            CPPUNIT_ASSERT(!hole.contains(Point(xout, yout)));
          }
      };

    const TriangulatorInterface::PolygonHole triangle(center, radius, 3);
    check_corners(triangle, 3);
    const TriangulatorInterface::PolygonHole diamond(center, radius, 4);
    check_corners(diamond, 4);
    const TriangulatorInterface::PolygonHole hexagon(center, radius, 6);
    check_corners(hexagon, 6);

    TriangulatorInterface::ArbitraryHole jaggy
    {{1,0}, {{0,-1},{2,-1},{2,1},{1.75,-.5},{1.5,1},{1.25,-.5},
             {1,1},{.75,-.5},{.5,1},{.25,-.5},{0,1}}};

    CPPUNIT_ASSERT(jaggy.contains({.1,-.3}));
    CPPUNIT_ASSERT(jaggy.contains({.5,.9}));
    CPPUNIT_ASSERT(jaggy.contains({.9,-.3}));
    CPPUNIT_ASSERT(jaggy.contains({1,0}));
    CPPUNIT_ASSERT(jaggy.contains({1.1,-.4}));
    CPPUNIT_ASSERT(jaggy.contains({1.5,.9}));
    CPPUNIT_ASSERT(jaggy.contains({1.9,-.3}));

    CPPUNIT_ASSERT(jaggy.contains({1.9,-.5}));
    CPPUNIT_ASSERT(jaggy.contains({1.6,-.5}));
    CPPUNIT_ASSERT(jaggy.contains({1.1,-.5}));
    CPPUNIT_ASSERT(jaggy.contains({.5,-.5}));
    CPPUNIT_ASSERT(jaggy.contains({.2,-.5}));

    CPPUNIT_ASSERT(!jaggy.contains({.1,.7}));
    CPPUNIT_ASSERT(!jaggy.contains({.5,1.1}));
    CPPUNIT_ASSERT(!jaggy.contains({.9,.7}));
    CPPUNIT_ASSERT(!jaggy.contains({1,-1.1}));
    CPPUNIT_ASSERT(!jaggy.contains({1.1,.8}));
    CPPUNIT_ASSERT(!jaggy.contains({1.5,1.1}));
    CPPUNIT_ASSERT(!jaggy.contains({1.9,.9}));

    CPPUNIT_ASSERT(!jaggy.contains({1.9,1}));
    CPPUNIT_ASSERT(!jaggy.contains({1.4,1}));
    CPPUNIT_ASSERT(!jaggy.contains({.9,1}));
    CPPUNIT_ASSERT(!jaggy.contains({.4,1}));
    CPPUNIT_ASSERT(!jaggy.contains({-.2,1}));
    CPPUNIT_ASSERT(!jaggy.contains({-.2,0}));
    CPPUNIT_ASSERT(!jaggy.contains({1.2,0}));

    TriangulatorInterface::ArbitraryHole square_jaggy
    {{1,0}, {{-.25,-1},{2,-1},{2,1},{1.75,1},{1.75,-.5},{1.5,-.5},
             {1.5,1},{1.25,1},{1.25,-.5},{1,-.5},{1,1},{.75,1},
             {.75,-.5},{.5,-.5},{.5,1},{.25,1},{.25,-.5},{0,-.5},
             {0,1},{-.25,1}}};

    CPPUNIT_ASSERT(square_jaggy.contains({-.1,-.3}));
    CPPUNIT_ASSERT(square_jaggy.contains({.4,.9}));
    CPPUNIT_ASSERT(square_jaggy.contains({.9,-.3}));
    CPPUNIT_ASSERT(square_jaggy.contains({.9,0}));
    CPPUNIT_ASSERT(square_jaggy.contains({1.1,-.6}));
    CPPUNIT_ASSERT(square_jaggy.contains({1.4,.9}));
    CPPUNIT_ASSERT(square_jaggy.contains({1.9,-.3}));

    CPPUNIT_ASSERT(!square_jaggy.contains({.1,-.3}));
    CPPUNIT_ASSERT(!square_jaggy.contains({.6,.9}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.1,-.3}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.1,0}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.1,-1.6}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.6,.9}));
    CPPUNIT_ASSERT(!square_jaggy.contains({2.1,-.3}));

    CPPUNIT_ASSERT(square_jaggy.contains({-.1,-.5}));
    CPPUNIT_ASSERT(square_jaggy.contains({.3,-.5}));
    CPPUNIT_ASSERT(square_jaggy.contains({.9,-.5}));
    CPPUNIT_ASSERT(square_jaggy.contains({1.3,-.5}));
    CPPUNIT_ASSERT(square_jaggy.contains({1.9,-.5}));

    CPPUNIT_ASSERT(!square_jaggy.contains({-.3,1}));
    CPPUNIT_ASSERT(!square_jaggy.contains({.2,1}));
    CPPUNIT_ASSERT(!square_jaggy.contains({.6,1}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.1,1}));
    CPPUNIT_ASSERT(!square_jaggy.contains({1.6,1}));
    CPPUNIT_ASSERT(!square_jaggy.contains({2.1,1}));
  }


  void testTriangulatorBase(MeshBase & mesh,
                            TriangulatorInterface & triangulator)
  {
    commonSettings(triangulator);

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), static_cast<dof_id_type>(2));
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


  void testTriangulatorTrapMesh (UnstructuredMesh & mesh)
  {
    // A non-square quad, so we don't have ambiguity about which
    // diagonal a Delaunay algorithm will pick.
    // Manually-numbered points, so we can use the point numbering as
    // a segment ordering even on DistributedMesh.
    mesh.add_point(Point(0,0), 0);
    mesh.add_point(Point(1,0), 1);
    mesh.add_point(Point(1,2), 2);
    mesh.add_point(Point(0,1), 3);
  }


  void testTriangulatorInterp(UnstructuredMesh & mesh,
                              TriangulatorInterface & triangulator,
                              int interpolate_boundary_points,
                              dof_id_type n_expected_elem,
                              Real expected_total_area = 1.5,
                              Real desired_area = 1000)
  {
    commonSettings(triangulator);

    if (!mesh.n_nodes())
      testTriangulatorTrapMesh(mesh);

    // Interpolate points!
    triangulator.set_interpolate_boundary_points(interpolate_boundary_points);

    // Try to insert points?
    triangulator.desired_area() = desired_area;

    triangulator.triangulate();

    if (n_expected_elem)
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

    commonSettings(triangulator);

    // Add a diamond hole in the center
    TriangulatorInterface::PolygonHole diamond(Point(0), std::sqrt(2)/2, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };
    triangulator.attach_hole_list(&holes);

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), static_cast<dof_id_type>(8));

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

    commonSettings(triangulator);

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

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), static_cast<dof_id_type>(12));

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


#ifdef LIBMESH_ENABLE_AMR
  void testTriangulatorRoundHole(MeshBase & mesh,
                                 TriangulatorInterface & triangulator)
  {
    // A square quad; we'll put a round hole in the middle.
    mesh.add_point(Point(19,19), 0);
    mesh.add_point(Point(21,19), 1);
    mesh.add_point(Point(21,21), 2);
    mesh.add_point(Point(19,21), 3);

    commonSettings(triangulator);

    // Add a square meshed hole in the center
    const Real radius = 0.5;
    const Point center{20,20};

    Mesh centermesh { mesh.comm() };
    MeshTools::Generation::build_sphere(centermesh, radius, 1, QUAD9, 0);
    MeshTools::Modification::translate(centermesh, center(0), center(1));

    TriangulatorInterface::MeshedHole centerhole { centermesh };

    CPPUNIT_ASSERT_EQUAL(centerhole.n_points(), 8u);
    Point inside = centerhole.inside();
    CPPUNIT_ASSERT_EQUAL(inside(0), Real(20));
    CPPUNIT_ASSERT_EQUAL(inside(1), Real(20));

    const std::vector<TriangulatorInterface::Hole*> holes { &centerhole };
    triangulator.attach_hole_list(&holes);
    triangulator.elem_type() = TRI6;

    triangulator.triangulate();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), static_cast<dof_id_type>(12));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(36));

    // Make sure we didn't screw up the outer sides.  We should
    // have exact values for the outer vertices, so we can use
    // those for a map.
    std::map<std::pair<Real, Real>, Point> outer_midpoints
      {{{19, 19}, {20, 19}},
       {{21, 19}, {21, 20}},
       {{21, 21}, {20, 21}},
       {{19, 21}, {19, 20}}
      };

    for (const auto & elem : mesh.active_local_element_ptr_range())
      for (const auto n : make_range(elem->n_sides()))
        {
          if (elem->neighbor_ptr(n))
            continue;

          auto it =
            outer_midpoints.find(std::make_pair(elem->point(n)(0),
                                                elem->point(n)(1)));
          if (it != outer_midpoints.end())
            {
              const Point error = it->second - elem->point(n+3);
              CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE,
                                  error.norm_sq());
            }
          else
            {
              const Point radius1 = elem->point(n) - center;
              CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE,
                                  std::abs(radius1.norm()-radius));

              const Point radius2 = elem->point((n+1)%3) - center;
              CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE,
                                  std::abs(radius2.norm()-radius));

              const Point radius3 = elem->point(n+3) - center;
              CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE,
                                  std::abs(radius3.norm()-radius));
            }
        }
  }
#endif // LIBMESH_ENABLE_AMR


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


  void testEdgesMesh(MeshBase & mesh, ElemType elem_type)
  {
    // The same quad as testTriangulator, but out of order
    auto node0 = mesh.add_point(Point(0,0), 0);
    auto node1 = mesh.add_point(Point(1,2), 1);
    auto node2 = mesh.add_point(Point(1,0), 2);
    auto node3 = mesh.add_point(Point(0,1), 3);

    // Edges, also out of order, but enough to put them in order
    auto edge13 = mesh.add_elem(Elem::build(elem_type));
    edge13->set_node(0) = node1;
    edge13->set_node(1) = node3;
    auto edge02 = mesh.add_elem(Elem::build(elem_type));
    edge02->set_node(0) = node0;
    edge02->set_node(1) = node2;
    auto edge30 = mesh.add_elem(Elem::build(elem_type));
    edge30->set_node(0) = node3;
    edge30->set_node(1) = node0;
    auto edge21 = mesh.add_elem(Elem::build(elem_type));
    edge21->set_node(0) = node2;
    edge21->set_node(1) = node1;

    // Add mid-edge nodes if asked to
    if (elem_type == EDGE3)
      {
        auto node4 = mesh.add_point(Point(.5,1.5), 4);
        edge13->set_node(2) = node4;
        auto node5 = mesh.add_point(Point(.5,0), 5);
        edge02->set_node(2) = node5;
        auto node6 = mesh.add_point(Point(0,.5), 6);
        edge30->set_node(2) = node6;
        auto node7 = mesh.add_point(Point(1,1), 7);
        edge21->set_node(2) = node7;
      }
    else
      libmesh_assert(elem_type == EDGE2);

    mesh.prepare_for_use();
  }


  void testHalfDomain(MeshBase & mesh,
                      TriangulatorInterface & triangulator,
                      ElemType elem_type)
  {
    // We might have 4 or 8 nodes on the outer boundary, depending on
    // whether it has mid-edge nodes
    const dof_id_type off = (elem_type == EDGE3)*4;

    // A pentagon we'll avoid via subdomain ids
    auto node4 = mesh.add_point(Point(2,0), 4+off);
    auto node5 = mesh.add_point(Point(3,0), 5+off);
    auto node6 = mesh.add_point(Point(3,2), 6+off);
    auto node7 = mesh.add_point(Point(2,2), 7+off);
    auto node8 = mesh.add_point(Point(2,1), 8+off);

    auto edge45 = mesh.add_elem(Elem::build(elem_type));
    edge45->set_node(0) = node4;
    edge45->set_node(1) = node5;
    edge45->subdomain_id() = 1;
    auto edge56 = mesh.add_elem(Elem::build(elem_type));
    edge56->set_node(0) = node5;
    edge56->set_node(1) = node6;
    edge56->subdomain_id() = 1;
    auto edge67 = mesh.add_elem(Elem::build(elem_type));
    edge67->set_node(0) = node6;
    edge67->set_node(1) = node7;
    edge67->subdomain_id() = 1;
    auto edge78 = mesh.add_elem(Elem::build(elem_type));
    edge78->set_node(0) = node7;
    edge78->set_node(1) = node8;
    edge78->subdomain_id() = 1;
    auto edge84 = mesh.add_elem(Elem::build(elem_type));
    edge84->set_node(0) = node8;
    edge84->set_node(1) = node4;
    edge84->subdomain_id() = 1;

    if (elem_type == EDGE3)
      {
        auto node9 = mesh.add_point(Point(2.5,0), 9+off);
        edge45->set_node(2) = node9;
        auto node10 = mesh.add_point(Point(3,1), 10+off);
        edge56->set_node(2) = node10;
        auto node11 = mesh.add_point(Point(2.5,2), 11+off);
        edge67->set_node(2) = node11;
        auto node12 = mesh.add_point(Point(2,1.5), 12+off);
        edge78->set_node(2) = node12;
        auto node13 = mesh.add_point(Point(2,.5), 13+off);
        edge84->set_node(2) = node13;
      }
    else
      libmesh_assert(elem_type == EDGE2);

    testEdgesMesh(mesh, elem_type);

    std::set<std::size_t> bdy_ids {0};
    triangulator.set_outer_boundary_ids(bdy_ids);

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


  void testTriangleHalfDomain()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testHalfDomain(mesh, triangle, EDGE2);
  }


  void testTriangleInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorInterp(mesh, triangle, 1, 6);
  }


  void testTriangleInterp2()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorInterp(mesh, triangle, 2, 10);
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


#ifdef LIBMESH_ENABLE_AMR
  void testTriangleRoundHole()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangle(mesh);
    testTriangulatorRoundHole(mesh, triangle);
  }
#endif // LIBMESH_ENABLE_AMR


  void testTriangleEdges()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    TriangleInterface triangulator(mesh);
    testEdgesMesh(mesh, EDGE2);

    this->testTriangulatorBase(mesh, triangulator);
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


  void testPoly2TriHalfDomain()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testHalfDomain(mesh, p2t_tri, EDGE2);
  }


  void testPoly2TriHalfDomainEdge3()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testHalfDomain(mesh, p2t_tri, EDGE3);
  }


  void testPoly2TriInterp()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorInterp(mesh, p2t_tri, 1, 6);
  }


  void testPoly2TriInterp2()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorInterp(mesh, p2t_tri, 2, 10);
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


#ifdef LIBMESH_ENABLE_AMR
  void testPoly2TriRoundHole()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorRoundHole(mesh, p2t_tri);
  }
#endif // LIBMESH_ENABLE_AMR


  void testPoly2TriEdges()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);
    testEdgesMesh(mesh, EDGE2);

    this->testTriangulatorBase(mesh, triangulator);
  }


  void testPoly2TriEdge3s()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);
    testEdgesMesh(mesh, EDGE3);

    this->testTriangulatorBase(mesh, triangulator);
  }


  void testPoly2TriBadEdges()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);

    // The same quad as testTriangulator, but out of order
    auto node0 = mesh.add_point(Point(0,0), 0);
    auto node1 = mesh.add_point(Point(1,2), 1);
    auto node2 = mesh.add_point(Point(1,0), 2);
    auto node3 = mesh.add_point(Point(0,1), 3);

    // Edges, but not enough to complete the quad
    auto edge13 = mesh.add_elem(Elem::build(EDGE2));
    edge13->set_node(0) = node1;
    edge13->set_node(1) = node3;
    auto edge02 = mesh.add_elem(Elem::build(EDGE2));
    edge02->set_node(0) = node0;
    edge02->set_node(1) = node2;
    auto edge30 = mesh.add_elem(Elem::build(EDGE2));
    edge30->set_node(0) = node3;
    edge30->set_node(1) = node0;

    mesh.prepare_for_use();

    testExceptionBase(mesh, triangulator, "Bad edge topology");
  }


  void testPoly2TriBad1DMultiBoundary()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);

    // Two separate triangles
    auto node0 = mesh.add_point(Point(0,0), 0);
    auto node1 = mesh.add_point(Point(0,1), 1);
    auto node2 = mesh.add_point(Point(1,0), 2);

    auto node3 = mesh.add_point(Point(2,0), 3);
    auto node4 = mesh.add_point(Point(2,1), 4);
    auto node5 = mesh.add_point(Point(3,0), 5);

    auto edge01 = mesh.add_elem(Elem::build(EDGE2));
    edge01->set_node(0) = node0;
    edge01->set_node(1) = node1;
    auto edge12 = mesh.add_elem(Elem::build(EDGE2));
    edge12->set_node(0) = node1;
    edge12->set_node(1) = node2;
    auto edge20 = mesh.add_elem(Elem::build(EDGE2));
    edge20->set_node(0) = node2;
    edge20->set_node(1) = node0;

    auto edge34 = mesh.add_elem(Elem::build(EDGE2));
    edge34->set_node(0) = node3;
    edge34->set_node(1) = node4;
    auto edge45 = mesh.add_elem(Elem::build(EDGE2));
    edge45->set_node(0) = node4;
    edge45->set_node(1) = node5;
    auto edge53 = mesh.add_elem(Elem::build(EDGE2));
    edge53->set_node(0) = node5;
    edge53->set_node(1) = node3;

    mesh.prepare_for_use();

    testExceptionBase(mesh, triangulator, "multiple loops of Edge");
  }

  void testPoly2TriBad2DMultiBoundary()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);

    // Two separate triangles
    auto node0 = mesh.add_point(Point(0,0), 0);
    auto node1 = mesh.add_point(Point(1,0), 1);
    auto node2 = mesh.add_point(Point(0,1), 2);

    auto node3 = mesh.add_point(Point(2,0), 3);
    auto node4 = mesh.add_point(Point(3,0), 4);
    auto node5 = mesh.add_point(Point(2,1), 5);

    auto tri012 = mesh.add_elem(Elem::build(TRI3));
    tri012->set_node(0) = node0;
    tri012->set_node(1) = node1;
    tri012->set_node(2) = node2;
    auto tri345 = mesh.add_elem(Elem::build(TRI3));
    tri345->set_node(0) = node3;
    tri345->set_node(1) = node4;
    tri345->set_node(2) = node5;

    mesh.prepare_for_use();

    testExceptionBase(mesh, triangulator, "cannot choose one");
  }


  void testPoly2TriEdgesRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator triangulator(mesh);
    testEdgesMesh(mesh, EDGE2);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 14);
  }


  void testPoly2TriSegments()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    Poly2TriTriangulator p2t_tri(mesh);
    testTriangulatorSegments(mesh, p2t_tri);
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

    commonSettings(triangulator);

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

  void testPoly2TriRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 15);
  }

  void testPoly2TriNonRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    // Make sure we see 0 as "don't refine", not "infinitely refine"
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 2, 0);
  }


  void testPoly2TriExtraRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 150, 0.01);
  }

  void testPoly2TriNonUniformRefined()
  {
    ParsedFunction<Real> var_area {"0.002*(1+2*x)*(1+2*y)"};
    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, nullptr, 1.5, 150, 0, &var_area);
  }

  void testPoly2TriHolesRefined()
  {
    LOG_UNIT_TEST;

    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 13);
  }

  void testPoly2TriHolesInterpRefined()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    Poly2TriTriangulator p2t_tri(mesh);

    Real total_area = 1.5;

    // Try a narrower point
    mesh.node_ref(2)(1) = 3;
    total_area += 0.5;

    // Add a bunch of tiny diamond holes
    const int N=3, M=3;
    TriangulatorInterface::PolygonHole diamond(Point(), std::sqrt(2)/20, 4);

    std::vector<TriangulatorInterface::AffineHole> hole_data;
    // Reserve so we don't invalidate pointers
    hole_data.reserve(M*N);
    std::vector<TriangulatorInterface::Hole*> holes;
    for (int i : make_range(M))
      for (int j : make_range(N))
        {
          Point shift(Real(i+1)/(M+1),Real(j+1)/(N+1));
          hole_data.emplace_back(diamond, 0, shift);
          holes.push_back(&hole_data.back());
          total_area -= hole_data.back().area();
        }

    p2t_tri.attach_hole_list(&holes);

    p2t_tri.set_refine_boundary_allowed(false);

    testTriangulatorInterp(mesh, p2t_tri, 4, 0, total_area, 0.03);
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
    testTriangulatorTrapMesh(mesh);
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
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 125, 0.01);
  }

  void testPoly2TriHolesNonUniformRefined()
  {
    // Add a diamond hole
    TriangulatorInterface::PolygonHole diamond(Point(0.5,0.5), std::sqrt(2)/4, 4);
    const std::vector<TriangulatorInterface::Hole*> holes { &diamond };

    ParsedFunction<Real> var_area {"0.002*(0.25+2*x)*(0.25+2*y)"};
    Mesh mesh(*TestCommWorld);
    testTriangulatorTrapMesh(mesh);
    testPoly2TriRefinementBase(mesh, &holes, 1.25, 150, 0, &var_area);
  }


#endif // LIBMESH_HAVE_POLY2TRI

};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshTriangulationTest );
