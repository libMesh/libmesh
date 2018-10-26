#include "test_comm.h"

#include <libmesh/bounding_box.h>
#include <tuple>
#include <algorithm>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

using namespace libMesh;

class BBoxTest : public CppUnit::TestCase {

public:
  CPPUNIT_TEST_SUITE( BBoxTest );
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( test_one_degenerate );
  CPPUNIT_TEST( test_two_degenerate );
  CPPUNIT_TEST( test_no_degenerate );
  CPPUNIT_TEST( test_signed_distance );
#endif
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void test_one_degenerate()
  {
    // Degenerate + Non-degenerate BBox tests: A unit square bounding
    // box is intersected with the following degenerate (planar)
    // bounding boxes:
    // 1.) well inside the unit bbox
    // 2.) just barely inside the unit bbox
    // 3.) on the surface of the unit bbox
    // 4.) just outside the surface of the unit bbox
    // 5.) well outside the unit bbox
    // For the "exact" intersection test: Tests 1-3 should return true, the rest false.
    // For the "fuzzy" intersection test: Tests 1-4 should return true with large enough TOLERANCE, the rest false.

    // The smallest number for which 1 and 1 + epsilon compare unequal.
    const Real eps = std::numeric_limits<Real>::epsilon();

    // Create unit square non-degenerate BBox for intersection testing.
    const Point
      min(0., 0., 0.),
      max(1., 1., 1.);
    BoundingBox non_degenerate(min, max);

    // Make (position, expected-exact-intersection, expected-fuzzy-intersection) tuples
    typedef std::tuple<Real, bool, bool> TestTuple;
    std::vector<TestTuple> tests =
      {
        std::make_tuple(0.5,       true,  true),  // inside
        std::make_tuple(1.0 - eps, true,  true),  // barely inside
        std::make_tuple(1.0,       true,  true),  // on
        std::make_tuple(1.0 + eps, false, true),  // barely outside
        std::make_tuple(1.5,       false, false)  // outside
      };

    for (unsigned int dir = 0; dir < 3; ++dir)
      for (const auto & t : tests)
        {
          // Create degenerate bounding box
          Point dmin = min, dmax = max;
          dmin(dir) = dmax(dir) = std::get<0>(t);
          BoundingBox degenerate(dmin, dmax);

          // std::cout << "degenerate.min() = " << degenerate.min() << std::endl;
          // std::cout << "degenerate.max() = " << degenerate.max() << std::endl;

          // Exact tests
          CPPUNIT_ASSERT(non_degenerate.intersects(degenerate) == std::get<1>(t));
          CPPUNIT_ASSERT(degenerate.intersects(non_degenerate) == std::get<1>(t));

          // Fuzzy tests
          CPPUNIT_ASSERT(non_degenerate.intersects(degenerate, /*abstol=*/TOLERANCE) == std::get<2>(t));
          CPPUNIT_ASSERT(degenerate.intersects(non_degenerate, /*abstol=*/TOLERANCE) == std::get<2>(t));
        }
  }

  void test_two_degenerate()
  {
    // Degenerate + Degenerate BBox tests: test intersections of unit
    // square degenerate (planar) bounding boxes in the x, y, and z
    // directions in the following cases:
    // 1.) Comparison BBox right on top of the original BBox.
    // 2.) Comparison BBox slightly to the "left" of the original BBox.
    // 3.) Comparison BBox slightly to the "right" of the original BBox.
    // 4.) Comparison BBox far away from the original BBox.

    // The smallest number for which 1 and 1 + epsilon compare unequal.
    const Real eps = std::numeric_limits<Real>::epsilon();
    // To test intersections of degenerate x, y, and z plane BBoxes,
    // we start with an x-plane and then permute the entries of
    // these vectors to subsequently test in the y and z directions.
    std::vector<Real>
      mins = {0.5, 0, 0},
      maxs = {.5, 1, 1};

    // Make (position, expected-exact-intersection, expected-fuzzy-intersection) tuples
    typedef std::tuple<Real, bool, bool> TestTuple;
    std::vector<TestTuple> tests =
      {
        std::make_tuple(0.5,       true,  true),  // on
        std::make_tuple(0.5 - eps, false, true),  // barely left
        std::make_tuple(0.5 + eps, false, true),  // barely right
        std::make_tuple(1.0,       false, false)  // off
      };

    for (unsigned int dir = 0; dir < 3; ++dir)
      {
        // Debugging
        // std::cout << "dir = " << dir
        //           << ", mins = " << mins[0] << ", " << mins[1] << ", " << mins[2]
        //           << ", maxs = " << maxs[0] << ", " << maxs[1] << ", " << maxs[2]
        //           << std::endl;

        const Point
          min(Point(mins[0], mins[1], mins[2])),
          max(Point(maxs[0], maxs[1], maxs[2]));
        BoundingBox initial(min, max);

        // Create each degenerate comparison BoundingBox and make
        // sure the expected intersection result is found.
        for (const auto & t : tests)
          {
            // Create comparison bbox
            Point cmin = min, cmax = max;
            cmin(dir) = cmax(dir) = std::get<0>(t);
            BoundingBox comparison(cmin, cmax);

            // Exact tests
            CPPUNIT_ASSERT(initial.intersects(comparison) == std::get<1>(t));
            CPPUNIT_ASSERT(comparison.intersects(initial) == std::get<1>(t));

            // Fuzzy tests
            CPPUNIT_ASSERT(initial.intersects(comparison, /*abstol=*/TOLERANCE) == std::get<2>(t));
            CPPUNIT_ASSERT(comparison.intersects(initial, /*abstol=*/TOLERANCE) == std::get<2>(t));
          }

        // Go to the next cyclic permutation of mins and maxs by
        // rotating the last entry off the end and onto the front of
        // the array.
        std::rotate(mins.rbegin(), mins.rbegin()+1, mins.rend());
        std::rotate(maxs.rbegin(), maxs.rbegin()+1, maxs.rend());
      }
  }

  void test_no_degenerate()
  {
    // The smallest number for which 1 and 1 + epsilon compare unequal.
    const Real eps = std::numeric_limits<Real>::epsilon();

    // Non-degenerate + Non-degenerate BBox tests. These tests
    // consider two initially unit-sized bounding boxes "stacked" in
    // the x, y, and z directions with slight perturbations to bring
    // them into intersection or not.
    std::vector<Real>
      mins = {1, 0, 0},
      maxs = {2, 1, 1};

    // All tests use the same unit bounding box for base comparisons.
    BoundingBox initial(Point(0.,0.,0.), Point(1.,1.,1.));

    // Make (position, expected-exact-intersection, expected-fuzzy-intersection) tuples
    typedef std::tuple<Real, bool, bool> TestTuple;
    std::vector<TestTuple> tests =
      {
        std::make_tuple(1.,       true,  true),  // on
        std::make_tuple(1. + eps, false, true),  // barely outside
        std::make_tuple(1. - eps, true,  true),  // barely inside
        std::make_tuple(1.5,      false, false), // definitely outside
        std::make_tuple(0.5,      true,  true)   // definitely inside
      };

    for (unsigned int dir = 0; dir < 3; ++dir)
      {
        // Debugging
        // std::cout << "dir = " << dir
        //           << ", mins = " << mins[0] << ", " << mins[1] << ", " << mins[2]
        //           << ", maxs = " << maxs[0] << ", " << maxs[1] << ", " << maxs[2]
        //           << std::endl;

        const Point
          min(Point(mins[0], mins[1], mins[2])),
          max(Point(maxs[0], maxs[1], maxs[2]));

        // Create each degenerate comparison BoundingBox and make
        // sure the expected intersection result is found.
        for (const auto & t : tests)
          {
            // Create comparison bbox with perturbed minimum coordinate.
            Point cmin = min, cmax = max;
            cmin(dir) = std::get<0>(t);
            BoundingBox comparison(cmin, cmax);

            // std::cout << "cmin = " << cmin << ", cmax = " << cmax << std::endl;

            // Exact tests
            CPPUNIT_ASSERT(initial.intersects(comparison) == std::get<1>(t));
            CPPUNIT_ASSERT(comparison.intersects(initial) == std::get<1>(t));

            // Fuzzy tests
            CPPUNIT_ASSERT(initial.intersects(comparison, /*abstol=*/TOLERANCE) == std::get<2>(t));
            CPPUNIT_ASSERT(comparison.intersects(initial, /*abstol=*/TOLERANCE) == std::get<2>(t));
          }

        // Go to the next cyclic permutation of mins and maxs by
        // rotating the last entry off the end and onto the front of
        // the array.
        std::rotate(mins.rbegin(), mins.rbegin()+1, mins.rend());
        std::rotate(maxs.rbegin(), maxs.rbegin()+1, maxs.rend());
      }
  }

  void test_signed_distance()
  {
    // A "unit" size bounding box for making distance comparisons.
    BoundingBox unit(Point(0.,0.,0.), Point(1.,1.,1.));

    // Test points inside the box
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, 0.5, 0.5)), -0.5, TOLERANCE * TOLERANCE);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, 0.6, 0.5)), -0.4, TOLERANCE * TOLERANCE);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.4, 0.5, 0.5)), -0.4, TOLERANCE * TOLERANCE);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.1, 0.1, 0.1)), -0.1, TOLERANCE * TOLERANCE);

    // Test points on the box
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(1.0, 0.5, 0.5)), 0., TOLERANCE * TOLERANCE);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(1.0, 0., 0.)), 0., TOLERANCE * TOLERANCE);

    // Test points outside the box
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(1.5, 0.5, 0.5)), 0.5, TOLERANCE * TOLERANCE);  // right
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(-0.5, 0.5, 0.5)), 0.5, TOLERANCE * TOLERANCE); // left
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, 0.5, 1.5)), 0.5, TOLERANCE * TOLERANCE);  // above
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, 0.5, -0.5)), 0.5, TOLERANCE * TOLERANCE); // below
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, -0.5, 0.5)), 0.5, TOLERANCE * TOLERANCE); // front
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(0.5, 1.5, 0.5)), 0.5, TOLERANCE * TOLERANCE);  // back

    // Outside the box, closest to a corner.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(2., 2., 2.)), std::sqrt(3.), TOLERANCE * TOLERANCE);    // Point along line (0,0,0) -> (1,1,1)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(-1., -1., -1.)), std::sqrt(3.), TOLERANCE * TOLERANCE); // Point along line (0,0,0) -> (1,1,1)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(1.5, 1.5, -0.5)), std::sqrt(3.)/2., TOLERANCE * TOLERANCE); // Point along line (0.5,0.5,0.5) -> (1,1,0)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(unit.signed_distance(Point(1.5, -0.5, -0.5)), std::sqrt(3.)/2., TOLERANCE * TOLERANCE); // Point along line (0.5,0.5,0.5) -> (1,0,0)
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( BBoxTest );
