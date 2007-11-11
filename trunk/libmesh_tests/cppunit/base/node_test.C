#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <node.h>

#include "../geom/point_test.h"

class NodeTest : public PointTestBase<Node> { 
public: 
  CPPUNIT_TEST_SUITE( NodeTest );

  POINTTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest );
