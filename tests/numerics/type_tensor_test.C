// libmesh includes
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>
#include <libmesh/point.h>

#include "libmesh_cppunit.h"


using namespace libMesh;

class TypeTensorTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  LIBMESH_CPPUNIT_TEST_SUITE(TypeTensorTest);

#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testInverse);
  CPPUNIT_TEST(testLeftMultiply);
  CPPUNIT_TEST(testRotation);
#endif
  CPPUNIT_TEST(testRowCol);
  CPPUNIT_TEST(testIsZero);
  CPPUNIT_TEST(testIsHPD);
#ifdef LIBMESH_HAVE_METAPHYSICL
  CPPUNIT_TEST(testReplaceAlgebraicType);
#endif
  CPPUNIT_TEST_SUITE_END();


private:
  void testInverse()
  {
    LOG_UNIT_TEST;

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
    LOG_UNIT_TEST;

    TensorValue<Real> tensor(1, 2, 0, 3, 4, 0);
    VectorValue<Real> vector(5, 6, 0);
    auto left_mult = vector * tensor;
    auto right_mult = tensor * vector;
    LIBMESH_ASSERT_FP_EQUAL(23, left_mult(0), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(34, left_mult(1), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(17, right_mult(0), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(39, right_mult(1), 1e-12);
  }

  void testOuterProduct()
  {
    LOG_UNIT_TEST;

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
    LOG_UNIT_TEST;

    {
      TensorValue<double> tensor;
      CPPUNIT_ASSERT(tensor.is_zero());
    }
    {
      TensorValue<double> tensor(0,1,2,3,4,5,6,7,8);
      CPPUNIT_ASSERT(!tensor.is_zero());
    }
  }

  void testIsHPD()
  {
    LOG_UNIT_TEST;

    {
#if LIBMESH_DIM == 3
      {
        TensorValue<double> tensor(2., 1., 0.,
                                   1., 2., 1.,
                                   0., 1., 2.);

        CPPUNIT_ASSERT(tensor.is_hpd(/*rel_tol=*/0.));
      }
      {
        // Symmetric but not positive-definite, because the upper 2x2
        // principal minor is zero.
        TensorValue<double> tensor(1,  0., 0.,
                                   0., 0., 1.,
                                   0., 1., 0.);

        CPPUNIT_ASSERT(!tensor.is_hpd());
      }
      {
        // Random matrix that is SPD by construction

        auto get_random = [&]()
        {
          return static_cast<Real>(std::rand()) / RAND_MAX; // in range [0,1]
        };

        TensorValue<double> tensor(get_random(), get_random(), get_random(),
                                   get_random(), get_random(), get_random(),
                                   get_random(), get_random(), get_random());

        // Make symmetric
        tensor = 0.5*(tensor + tensor.transpose());

        // Make positive definite
        for (auto i : make_range(LIBMESH_DIM))
          tensor(i,i) += 3.;

        CPPUNIT_ASSERT(tensor.is_hpd());
      }

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

      auto i = Complex(0, 1);
      auto one = Complex(1, 0);
      {
        // Symmetric but not Hermitian, because Hermitian matrices must
        // have real-valued entries on the diagonal.
        TensorValue<Complex> tensor(i,  0., 0.,
                                    0., 1., 0.,
                                    0., 0., 1.);

        CPPUNIT_ASSERT(!tensor.is_hpd());
      }
      {
        // Symmetric but not Hermitian, because the off-diagonal
        // entries must be complex conjugates of one another.
        TensorValue<Complex> tensor(2,  i,  0.,
                                    i,  2., i,
                                    0., i,  2.);

        CPPUNIT_ASSERT(!tensor.is_hpd());
      }
      {
        // Both Hermitian and positive-definite
        TensorValue<Complex> tensor(2.,     one+i,  0.,
                                    one-i,  2.,     one+i,
                                    0.,     one-i,  2.);

        CPPUNIT_ASSERT(tensor.is_hpd());
      }

#endif

#endif // LIBMESH_DIM == 3
    }
  }

  void testRotation()
  {
    LOG_UNIT_TEST;

    {
      Point x(1, 0, 0);
      const auto R = RealTensorValue::extrinsic_rotation_matrix(90, 0, 0);
      auto rotated = R * x;
      constexpr auto tol = TOLERANCE * TOLERANCE;
      LIBMESH_ASSERT_FP_EQUAL(0, rotated(0), tol);
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(1), tol);
      LIBMESH_ASSERT_FP_EQUAL(0, rotated(2), tol);

      const auto invR = RealTensorValue::inverse_extrinsic_rotation_matrix(90, 0, 0);
      rotated = invR * rotated;
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(0), tol);
      LIBMESH_ASSERT_FP_EQUAL(0, rotated(1), tol);
      LIBMESH_ASSERT_FP_EQUAL(0, rotated(2), tol);
    }

    {
      Point x(1, 1, 1);
      const auto R = RealTensorValue::extrinsic_rotation_matrix(90, 90, 90);
      auto rotated = R * x;

      constexpr auto tol = TOLERANCE * TOLERANCE;
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(0), tol);
      LIBMESH_ASSERT_FP_EQUAL(-1, rotated(1), tol);
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(2), tol);

      const auto invR = RealTensorValue::inverse_extrinsic_rotation_matrix(90, 90, 90);
      rotated = invR * rotated;
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(0), tol);
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(1), tol);
      LIBMESH_ASSERT_FP_EQUAL(1, rotated(2), tol);
    }
  }

  void testRowCol()
  {
    LOG_UNIT_TEST;

    TensorValue<Real> t;
    for (unsigned int i = 0; i < LIBMESH_DIM; i++)
      for (unsigned int j = 0; j < LIBMESH_DIM; j++)
        t(i, j) = Real(i * LIBMESH_DIM + j);

    constexpr Real tol = TOLERANCE * TOLERANCE;
    for (unsigned int k = 0; k < LIBMESH_DIM; k++)
    {
      const auto row = t.row(k);
      for (unsigned int l = 0; l < LIBMESH_DIM; l++)
        LIBMESH_ASSERT_FP_EQUAL(Real(k * LIBMESH_DIM + l), row(l), tol);

      const auto col = t.column(k);
      for (unsigned int l = 0; l < LIBMESH_DIM; l++)
        LIBMESH_ASSERT_FP_EQUAL(Real(l * LIBMESH_DIM + k), col(l), tol);
    }
  }

#ifdef LIBMESH_HAVE_METAPHYSICL
  void testReplaceAlgebraicType()
  {
    typedef typename MetaPhysicL::ReplaceAlgebraicType<
        std::vector<TypeTensor<double>>,
        typename TensorTools::IncrementRank<
            typename MetaPhysicL::ValueType<std::vector<TypeTensor<double>>>::type>::type>::type
        ReplacedType;
    constexpr bool assertion =
        std::is_same<ReplacedType, std::vector<TypeNTensor<3,double>>>::value;
    CPPUNIT_ASSERT(assertion);
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION(TypeTensorTest);
