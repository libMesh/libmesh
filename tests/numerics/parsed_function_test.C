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
  CPPUNIT_TEST(testInlineGetter);
  CPPUNIT_TEST(testInlineSetter);

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

  void testInlineGetter()
  {
    ParsedFunction<Number> ax2("a:=4.5;a*x*2");

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(ax2(Point(0.25,0.25,0.25))), 2.25, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(ax2.get_inline_value("a")), 4.5, TOLERANCE*TOLERANCE);

    ParsedFunction<Number> cxy8
      ("a := 4 ; b := a/2+1; c:=b-a+3.5; c*x*2*y*4");

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8(Point(0.5,0.5,0.5))), 5.0, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8.get_inline_value("b")), 3.0, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8.get_inline_value("c")), 2.5, TOLERANCE*TOLERANCE);
  }

  void testInlineSetter()
  {
    ParsedFunction<Number> ax2("a:=4.5;a*x*2");
    ax2.set_inline_value("a", 2.5);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(ax2(Point(0.25,0.25,0.25))), 1.25, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(ax2.get_inline_value("a")), 2.5, TOLERANCE*TOLERANCE);

    ParsedFunction<Number> cxy8
      ("a := 4 ; b := a/2+1; c:=b-a+3.5; c*x*2*y*4");
    cxy8.set_inline_value("a", 2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8(Point(0.5,0.5,0.5))), 7.0, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8.get_inline_value("b")), 2.0, TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL
      (libmesh_real(cxy8.get_inline_value("c")), 3.5, TOLERANCE*TOLERANCE);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(ParsedFunctionTest);
