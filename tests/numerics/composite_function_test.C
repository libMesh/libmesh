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
#include <libmesh/analytic_function.h>

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
#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testTimeDependence);
#endif

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

      // Test that ConstFunction copy construction works.
      ConstFunction<Real> cf_one(1);
      ConstFunction<Real> cf_one_copy(cf_one);
      composite_inner.attach_subfunction (cf_one_copy, index_sets[0]);

      // Test that ConstFunction move construction works.
      ConstFunction<Real> cf_two(2);
      ConstFunction<Real> cf_two_move(std::move(cf_two));
      composite_inner.attach_subfunction (cf_two_move, index_sets[1]);

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
    // Test that ConstFunction copy- and move-assignment works.
    ConstFunction<Real> cf_three(3);
    ConstFunction<Real> cf_three_copy_assign(0);
    ConstFunction<Real> cf_three_move_assign(0);
    cf_three_copy_assign = cf_three;
    cf_three_move_assign = std::move(cf_three_copy_assign);
    composite_outer.attach_subfunction(cf_three_move_assign, index_sets[2]);

    // Test that move ctor works. Note that composite_outer should not
    // be used for anything once it has been moved from!
    CompositeFunction<Real> composite_outer_copy1(std::move(composite_outer));

    // Test that move assignment also works. The first copy should not be
    // used again after being move assigned.
    CompositeFunction<Real> composite_outer_copy2;
    composite_outer_copy2 = std::move(composite_outer_copy1);

    DenseVector<Real> test_two(8);
    composite_outer_copy2(Point(0), 0, test_two);

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

    // Test constructing AnalyticFunction with lambda.
    auto af_lambda =
      [](const Point & p, const Real t) -> Real
      { return p(0)*p(0) + p(1)*p(1) + t*t; };
    AnalyticFunction<Real> x2y2t2(af_lambda);

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

      // Test AnalyticFunction copy ctor and copy assignment
      AnalyticFunction<Real> x2y2t2_copy1(x2y2t2);
      AnalyticFunction<Real> x2y2t2_copy2([](const Point &, const Real) -> Real { return 0; });
      x2y2t2_copy2 = x2y2t2_copy1;
      composite.attach_subfunction(x2y2t2_copy2, index_set);
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

      // Test AnalyticFunction move ctor and move assignment. Note: we
      // first copy and then steal the copy's resources to avoid
      // messing with any later tests of "x2y2t2".
      AnalyticFunction<Real> x2y2t2_copy(x2y2t2);
      AnalyticFunction<Real> x2y2t2_move1(std::move(x2y2t2_copy));
      AnalyticFunction<Real> x2y2t2_move2([](const Point &, const Real) -> Real { return 0; });
      x2y2t2_move2 = std::move(x2y2t2_move1);
      composite.attach_subfunction(x2y2t2_move2, index_set);
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
