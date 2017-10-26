// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include <libmesh/composite_function.h>
#include <libmesh/const_function.h>
#include <libmesh/dense_vector.h>
#include <libmesh/parsed_function.h>
#include <libmesh/zero_function.h>

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

class CompositeFunctionTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(CompositeFunctionTest);

  CPPUNIT_TEST(testRemap);
  CPPUNIT_TEST(testTimeDependence);

  CPPUNIT_TEST_SUITE_END();


private:
  void testRemap()
  {
    std::vector<std::vector<unsigned int>> index_sets(4);
    index_sets[0].resize(2);
    index_sets[0][0] = 3;
    index_sets[0][1] = 4;
    index_sets[1].resize(3);
    index_sets[1][0] = 0;
    index_sets[1][1] = 1;
    index_sets[1][2] = 2;
    index_sets[2].resize(3);
    index_sets[2][0] = 0;
    index_sets[2][1] = 2;
    index_sets[2][2] = 4;
    index_sets[3].resize(5);
    index_sets[3][0] = 5;
    index_sets[3][1] = 1;
    index_sets[3][2] = 3;
    index_sets[3][3] = 6;
    index_sets[3][4] = 7;

    CompositeFunction<Real> composite_outer;

    {
      CompositeFunction<Real> composite_inner;
      composite_inner.attach_subfunction
        (ConstFunction<Real>(1), index_sets[0]);
      composite_inner.attach_subfunction
        (ConstFunction<Real>(2), index_sets[1]);
      composite_outer.attach_subfunction
        (composite_inner, index_sets[3]);

      DenseVector<Real> test_one(5);

      composite_inner(Point(0), 0, test_one);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(test_one(0), 2, 1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(test_one(1), 2, 1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(test_one(2), 2, 1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(test_one(3), 1, 1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(test_one(4), 1, 1.e-12);
    }
    composite_outer.attach_subfunction
      (ConstFunction<Real>(3), index_sets[2]);

    DenseVector<Real> test_two(8);
    composite_outer(Point(0), 0, test_two);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(0), 3, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(2), 3, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(4), 3, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(5), 2, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(1), 2, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(3), 2, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(6), 1, 1.e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_two(7), 1, 1.e-12);
  }

  void testTimeDependence()
  {

    // We'll test the order of adding these functions to
    // make sure time dependence gets detected/updated correctly
    // for each
    ParsedFunction<Real> no_t("x*2+y^2-tanh(z)+atan(x-y)");
    ParsedFunction<Real> no_t2("x*2+y^2+z^2");
    ZeroFunction<Real> zero;

    ParsedFunction<Real> xyt("x+y+t");
    ParsedFunction<Real> x2y2t2("x*2+y^2+t^2");


    std::vector<unsigned int> index_set(1,0);

    {
      // composite should not be time dependent since this is the first subfunction
      // added and it's not time-dependent
      CompositeFunction<Real> composite;
      composite.attach_subfunction(no_t, index_set);
      CPPUNIT_ASSERT(!composite.is_time_dependent());

      // Now composite should be time-dependent since we've now added a time dependent function
      index_set[0] = 1;
      composite.attach_subfunction(xyt, index_set);
      CPPUNIT_ASSERT(composite.is_time_dependent());

      // Composite should still be time-dependent
      index_set[0] = 2;
      composite.attach_subfunction(x2y2t2, index_set);
      CPPUNIT_ASSERT(composite.is_time_dependent());
    }


    {
      CompositeFunction<Real> composite;

      // composite should be time-dependent since we've added a time dependent function
      index_set[0] = 0;
      composite.attach_subfunction(xyt, index_set);
      CPPUNIT_ASSERT(composite.is_time_dependent());

      // composite should still be time-dependent since the previous function was time-dependent
      index_set[0] = 1;
      composite.attach_subfunction(no_t, index_set);
      CPPUNIT_ASSERT(composite.is_time_dependent());

      // Composite should still be time-dependent
      index_set[0] = 2;
      composite.attach_subfunction(x2y2t2, index_set);
      CPPUNIT_ASSERT(composite.is_time_dependent());
    }

    {
      CompositeFunction<Real> composite;

      // composite should not be time-dependent since we've added a time independent function
      index_set[0] = 0;
      composite.attach_subfunction(no_t, index_set);
      CPPUNIT_ASSERT(!composite.is_time_dependent());

      // composite should still be time-independent
      index_set[0] = 1;
      composite.attach_subfunction(no_t2, index_set);
      CPPUNIT_ASSERT(!composite.is_time_dependent());

      // Composite should still be time-independent
      index_set[0] = 2;
      composite.attach_subfunction(zero, index_set);
      CPPUNIT_ASSERT(!composite.is_time_dependent());
    }

  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(CompositeFunctionTest);
