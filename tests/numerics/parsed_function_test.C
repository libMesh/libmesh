// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_function.h"
#include "libmesh/system.h"

#ifdef LIBMESH_HAVE_FPARSER

// test includes
#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

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
  CPPUNIT_TEST(testTimeDependence);

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

  void testTimeDependence()
  {
    ParsedFunction<Number> no_t("x*2+y^2-tanh(z)+atan(x-y)");
    CPPUNIT_ASSERT(!no_t.is_time_dependent());

    ParsedFunction<Number> xyt("x+y+t");
    CPPUNIT_ASSERT(xyt.is_time_dependent());

    ParsedFunction<Number> x2y2t2("x*2+y^2+t^2");
    CPPUNIT_ASSERT(x2y2t2.is_time_dependent());

    ParsedFunction<Number> ztanht("z*tanh(t)");
    CPPUNIT_ASSERT(ztanht.is_time_dependent());
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(ParsedFunctionTest);

#endif // #ifdef LIBMESH_HAVE_FPARSER
