// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

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

class TypeTensorTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(TypeTensorTest);

#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testInverse);
  CPPUNIT_TEST(testLeftMultiply);
#endif
  CPPUNIT_TEST(testIsZero);

  CPPUNIT_TEST_SUITE_END();


private:
  void testInverse()
  {
    // This random input tensor and its inverse came from Octave/Matlab:
    // > format long e
    // > A = rand(3)
    // > inv(A)

    // The base class, TypeTensor, has a protected constructor.  We
    // are using the derived class, TensorValue, for our tests...
    TensorValue<double> tensor(9.08973348886179e-01, 3.36455579239923e-01, 5.16389236893863e-01,
                               9.44156071777472e-01, 1.35610910092516e-01, 1.49881119060538e-02,
                               1.15988384086146e-01, 6.79845197685518e-03, 3.77028969454745e-01);

    TensorValue<double> inverse = tensor.inverse();

    TensorValue<double> true_inverse(-6.57484735104482e-01, 1.58926633961497e+00,  8.37330721137561e-01,
                                     4.56430940967411e+00, -3.64404559823061e+00, -6.10654107858520e+00,
                                     1.19965194510943e-01, -4.23210359257434e-01,  2.50483242797707e+00);

    for (unsigned i=0; i<3; ++i)
      for (unsigned j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(inverse(i,j), true_inverse(i,j), TOLERANCE * TOLERANCE);
  }

  void testLeftMultiply()
  {
    TensorValue<Real> tensor(1, 2, 0, 3, 4, 0);
    VectorValue<Real> vector(5, 6, 0);
    auto left_mult = vector * tensor;
    auto right_mult = tensor * vector;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(left_mult(0), 23, 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(left_mult(1), 34, 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(right_mult(0), 17, 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(right_mult(1), 39, 1e-12);
  }

  void testOuterProduct()
  {
    auto tol = TOLERANCE * TOLERANCE;
    VectorValue<Real> a(2, 3, 4);
    VectorValue<Real> b(5, 6, 7);
    auto product = outer_product(a, b);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(0, 0), 10, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(0, 1), 12, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(0, 2), 14, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(1, 0), 15, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(1, 1), 18, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(1, 2), 21, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(2, 0), 20, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(2, 1), 24, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(product(2, 2), 28, tol);
  }

  void testIsZero()
  {
    {
      TensorValue<double> tensor;
      CPPUNIT_ASSERT(tensor.is_zero());
    }
    {
      TensorValue<double> tensor(0,1,2,3,4,5,6,7,8);
      CPPUNIT_ASSERT(!tensor.is_zero());
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TypeTensorTest);
