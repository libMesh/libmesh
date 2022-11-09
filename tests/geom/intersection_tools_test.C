#include "libmesh_cppunit.h"

#include <libmesh/intersection_tools.h>
#include <libmesh/point.h>
#include <libmesh/int_range.h>

using namespace libMesh;

class IntersectionToolsTest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( IntersectionToolsTest );
  CPPUNIT_TEST( within_segment );
  CPPUNIT_TEST( collinear );
  CPPUNIT_TEST_SUITE_END();

public:

  void within_segment()
  {
    LOG_UNIT_TEST;

    const Point s1(1.0, 2.0, 3.0);
    const Point s2(2.0, 3.0, 4.0);
    const auto length_vec = s2 - s1;
    const auto length = length_vec.norm();
    const auto s1_to_s2 = length_vec / length;

    int segments = 10;
    Real dx = (Real)1 / segments * length;
    for (const auto i : make_range(-1, segments + 1))
    {
      const auto p = s1 + Real(i) * dx * s1_to_s2;
      IntersectionTools::WithinSegmentResult within_result = IntersectionTools::WithinSegmentResult::NOT_WITHIN;
      if (i == 0)
        within_result = IntersectionTools::WithinSegmentResult::AT_BEGINNING;
      else if (i > 0 && i < segments)
        within_result = IntersectionTools::WithinSegmentResult::BETWEEN;
      else if (i == segments)
        within_result = IntersectionTools::WithinSegmentResult::AT_END;

      CPPUNIT_ASSERT_EQUAL(IntersectionTools::within_segment(s1, s2, length, p), within_result);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::within_segment(s1, s2, p), within_result);
    }

    CPPUNIT_ASSERT_EQUAL(IntersectionTools::within_segment(s1, s2, Point(9.9, 5, 3)),
                         IntersectionTools::WithinSegmentResult::NOT_WITHIN);
  }

  void collinear()
  {
    LOG_UNIT_TEST;

    const auto do_assert = [](const Point & p1,
                              const Point & p2,
                              const Point & p3,
                              const bool collinear)
    {
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p1, p2, p3), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p1, p3, p2), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p2, p1, p3), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p2, p3, p1), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p3, p1, p2), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p3, p2, p1), collinear);
    };

    // two of the same points
    do_assert(Point(1, 2, 3), Point(1, 2, 3), Point(4, 5, 6), true);

    // three of the same points
    do_assert(Point(1, 2, 3), Point(1, 2, 3), Point(1, 2, 3), true);

    // all in a line
    do_assert(Point(1, 1, 2), Point(1, 1, 3), Point(1, 1, 0), true);

    // not collinear
    do_assert(Point(0, 1, 2), Point(0, 1, 3), Point(1, 5, 10), false);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( IntersectionToolsTest );
