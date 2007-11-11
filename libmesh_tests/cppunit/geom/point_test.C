#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <point.h>

#include "../numerics/type_vector_test.h"

class PointTest : public TypeVectorTestBase<Point> { 
public: 
  CPPUNIT_TEST_SUITE( PointTest );

  TYPEVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( PointTest );
