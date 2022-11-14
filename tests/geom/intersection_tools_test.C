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

};

CPPUNIT_TEST_SUITE_REGISTRATION( IntersectionToolsTest );
