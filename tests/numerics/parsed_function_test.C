// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_function.h"
#include "libmesh/system.h"

// test includes
#include "test_comm.h"

using namespace libMesh;

class ParsedFunctionTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(ParsedFunctionTest);

  CPPUNIT_TEST(testValues);

  CPPUNIT_TEST_SUITE_END();


private:

  void testValues()
  {
    ParsedFunction<Number> x2("x*2");

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(x2(Point(0.5,1.5,2.5))), 1.0, TOLERANCE*TOLERANCE);

    ParsedFunction<Number> xy8("x*y*8");

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(xy8(Point(0.5,1.5,2.5))), 6.0, TOLERANCE*TOLERANCE);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION(ParsedFunctionTest);
