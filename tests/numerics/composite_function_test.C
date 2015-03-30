// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include <libmesh/composite_function.h>
#include <libmesh/const_function.h>
#include <libmesh/dense_vector.h>

using namespace libMesh;

class CompositeFunctionTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(CompositeFunctionTest);

  CPPUNIT_TEST(testRemap);

  CPPUNIT_TEST_SUITE_END();


private:
  void testRemap()
  {
    std::vector<std::vector<unsigned int> > index_sets(4);
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

};

CPPUNIT_TEST_SUITE_REGISTRATION(CompositeFunctionTest);
