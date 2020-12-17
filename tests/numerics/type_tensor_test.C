// libmesh includes
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

// C++ includes
#include <cstdlib> // std::rand

#include "libmesh_cppunit.h"

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
  CPPUNIT_TEST(matMult3);
  CPPUNIT_TEST(axpy);
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
        LIBMESH_ASSERT_FP_EQUAL(inverse(i,j), true_inverse(i,j), 1e-12);
  }

  void testLeftMultiply()
  {
    TensorValue<Real> tensor(1, 2, 0, 3, 4, 0);
    VectorValue<Real> vector(5, 6, 0);
    auto left_mult = vector * tensor;
    auto right_mult = tensor * vector;
    LIBMESH_ASSERT_FP_EQUAL(23, left_mult(0), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(34, left_mult(1), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(17, right_mult(0), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(39, right_mult(1), 1e-12);
  }

  void matMult3()
  {
    // Test that the operations
    // .) A * B * C;
    // .) A *= B; A *= C;
    // .) mat_mult3(A, B, C);
    // all give the same result.
    auto r = [](){ return static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX); };

    RealTensorValue A(r(), r(), r(), r(), r(), r(), r(), r(), r());
    RealTensorValue B(r(), r(), r(), r(), r(), r(), r(), r(), r());
    RealTensorValue C(r(), r(), r(), r(), r(), r(), r(), r(), r());

    RealTensorValue D1 = A * B * C;
    RealTensorValue D2 = A; D2 *= B; D2 *= C;
    RealTensorValue D3 = RealTensorValue::mat_mult3(A, B, C);

    Real D12 = (D1 - D2).norm();
    Real D23 = (D2 - D3).norm();

    LIBMESH_ASSERT_FP_EQUAL(0, D12, TOLERANCE * TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(0, D23, TOLERANCE * TOLERANCE);
  }

  void axpy()
  {
    // Test that the operations
    // .) A * x + y;
    // .) axpy(A, x, y);
    // all give the same result.
    auto r = [](){ return static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX); };

    RealTensorValue A(r(), r(), r(), r(), r(), r(), r(), r(), r());
    RealVectorValue x(r(), r(), r());
    RealVectorValue y(r(), r(), r());

    RealVectorValue z1 = A * x + y;
    RealVectorValue z2 = RealTensorValue::axpy(A, x, y);

    Real z12 = (z1 - z2).norm();
    LIBMESH_ASSERT_FP_EQUAL(0, z12, TOLERANCE * TOLERANCE);
  }

  void testOuterProduct()
  {
    auto tol = TOLERANCE * TOLERANCE;
    VectorValue<Real> a(2, 3, 4);
    VectorValue<Real> b(5, 6, 7);
    auto product = outer_product(a, b);
    LIBMESH_ASSERT_FP_EQUAL(10, product(0, 0), tol);
    LIBMESH_ASSERT_FP_EQUAL(12, product(0, 1), tol);
    LIBMESH_ASSERT_FP_EQUAL(14, product(0, 2), tol);
    LIBMESH_ASSERT_FP_EQUAL(15, product(1, 0), tol);
    LIBMESH_ASSERT_FP_EQUAL(18, product(1, 1), tol);
    LIBMESH_ASSERT_FP_EQUAL(21, product(1, 2), tol);
    LIBMESH_ASSERT_FP_EQUAL(20, product(2, 0), tol);
    LIBMESH_ASSERT_FP_EQUAL(24, product(2, 1), tol);
    LIBMESH_ASSERT_FP_EQUAL(28, product(2, 2), tol);
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
