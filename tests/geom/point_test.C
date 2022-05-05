#include <libmesh/point.h>

#include "point_test.h"

using namespace libMesh;

class PointTest : public PointTestBase<Point> {
public:
  PointTest() :
    PointTestBase<Point>() {
    this->libmesh_suite_name = "PointTest";
  }

  CPPUNIT_TEST_SUITE( PointTest );

  POINTTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( PointTest );
