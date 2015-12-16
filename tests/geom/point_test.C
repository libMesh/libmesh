// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/point.h>

#include "point_test.h"

using namespace libMesh;

class PointTest : public PointTestBase<Point> {
public:
  CPPUNIT_TEST_SUITE( PointTest );

  POINTTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( PointTest );
