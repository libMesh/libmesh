// libMesh includes
#include "libmesh/rb_parameters.h"

// CPPUnit includes
#include "libmesh_cppunit.h"

using namespace libMesh;

class RBParametersTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE ( RBParametersTest );
  CPPUNIT_TEST( testScalar );
  CPPUNIT_TEST_SUITE_END();

public:

  // virtual void setUp()
  // {}

  // virtual void tearDown()
  // {}

  void testScalar()
  {
    LOG_UNIT_TEST;

    // Test scalar-valued interfaces
    RBParameters params;
    params.set_value("a", 1.);
    params.set_value("b", 2.);
    params.set_value("c", 3.);

    // Expected result
    // a: 1.000000e+00
    // b: 2.000000e+00
    // c: 3.000000e+00
    CPPUNIT_ASSERT(params.has_value("a"));
    CPPUNIT_ASSERT(params.has_value("b"));
    CPPUNIT_ASSERT(params.has_value("c"));
    CPPUNIT_ASSERT_EQUAL(params.get_value("a"), 1.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("b"), 2.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("c"), 3.);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION ( RBParametersTest );
