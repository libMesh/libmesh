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
#include "libmesh_cppunit.h"


using namespace libMesh;

class ParsedFunctionTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  LIBMESH_CPPUNIT_TEST_SUITE(ParsedFunctionTest);

#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testValues);
  CPPUNIT_TEST(testInlineGetter);
  CPPUNIT_TEST(testInlineSetter);
  CPPUNIT_TEST(testTimeDependence);
#endif

  CPPUNIT_TEST_SUITE_END();


private:

  void testValues()
  {
    LOG_UNIT_TEST;

    // Test that we can copy these into vectors
    std::vector<ParsedFunction<Number>> pfvec;

    {
      ParsedFunction<Number> xy8("x*y*8");

      // Test that the move ctor works
      ParsedFunction<Number> xy8_stolen(std::move(xy8));

      pfvec.push_back(xy8_stolen);

      LIBMESH_ASSERT_NUMBERS_EQUAL
        (6.0, xy8_stolen(Point(0.5,1.5,2.5)), TOLERANCE*TOLERANCE);
    }
    LIBMESH_ASSERT_NUMBERS_EQUAL
      (6.0, pfvec[0](Point(0.5,1.5,2.5)), TOLERANCE*TOLERANCE);
  }

  void testInlineGetter()
  {
    LOG_UNIT_TEST;

    ParsedFunction<Number> ax2("a:=4.5;a*x*2");

    // Test whether move assignment works.
    ParsedFunction<Number> ax2_stolen("x");
    ax2_stolen = std::move(ax2);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (2.25, ax2_stolen(Point(0.25,0.25,0.25)), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (4.5, ax2_stolen.get_inline_value("a"), TOLERANCE*TOLERANCE);

    ParsedFunction<Number> cxy8
      ("a := 4 ; b := a/2+1; c:=b-a+3.5; c*x*2*y*4");

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (5.0, cxy8(Point(0.5,0.5,0.5)), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (3.0, cxy8.get_inline_value("b"), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (2.5, cxy8.get_inline_value("c"), TOLERANCE*TOLERANCE);
  }

  void testInlineSetter()
  {
    LOG_UNIT_TEST;

    ParsedFunction<Number> ax2("a:=4.5;a*x*2");
    ax2.set_inline_value("a", 2.5);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (1.25, ax2(Point(0.25,0.25,0.25)), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (2.5, ax2.get_inline_value("a"), TOLERANCE*TOLERANCE);

    ParsedFunction<Number> cxy8
      ("a := 4 ; b := a/2+1; c:=b-a+3.5; c*x*2*y*4");
    cxy8.set_inline_value("a", 2);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (7.0, cxy8(Point(0.5,0.5,0.5)), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (2.0, cxy8.get_inline_value("b"), TOLERANCE*TOLERANCE);

    LIBMESH_ASSERT_NUMBERS_EQUAL
      (3.5, cxy8.get_inline_value("c"), TOLERANCE*TOLERANCE);
  }

  void testTimeDependence()
  {
    LOG_UNIT_TEST;

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
