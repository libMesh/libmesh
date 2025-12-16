// libmesh includes
#include <libmesh/type_n_tensor.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>
#include <libmesh/point.h>

#include "libmesh_cppunit.h"

using namespace libMesh;

class TypeNTensorTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  LIBMESH_CPPUNIT_TEST_SUITE(TypeNTensorTest);

#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testOperatorsScalar);
  CPPUNIT_TEST(testOperatorsTensor);
  CPPUNIT_TEST(testCastVector);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testNorm);
  CPPUNIT_TEST(testSlice);
#endif

  CPPUNIT_TEST_SUITE_END();

private:
  void testOperatorsScalar()
  {
    LOG_UNIT_TEST;

    std::vector<Real> values = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    TypeNTensor<2, Real> T(values);

    // Multiply by a number
    T *= 3.1;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(3.1 * values[i], T._coords[i], 1e-12);

    // Divide by a number
    T /= 3.1;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i], T._coords[i], 1e-12);
  }

  void testOperatorsTensor()
  {
    std::vector<Real> values = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Real> values2 = {0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
    TypeNTensor<2, Real> T(values);
    TypeNTensor<2, Real> T2(values2);

    // Add a tensor
    T += T2;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i] + values2[i], T._coords[i], 1e-12);

    // Return sum of two tensors
    const auto T3 = T + T2;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i] + 2 * values2[i], T3._coords[i], 1e-12);

    // Subtract a tensor
    T -= T2;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i], T._coords[i], 1e-12);

    // Return subtraction of one tensor by another
    const auto T4 = T - T2;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i] - values2[i], T4._coords[i], 1e-12);

    // Add a scaled tensor
    T.add_scaled(T2, 3);
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i] + 3 * values2[i], T._coords[i], 1e-12);
  }

  void testCastVector()
  {
    std::vector<Real> values = {0, 1, 2};
    TypeNTensor<1, Real> T(values);
    VectorValue<Real> v2 = T;

    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(values[i], v2(i), 1e-12);
  }

  void testZero()
  {
    std::vector<Real> values = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    TypeNTensor<2, Real> T(values);
    TypeNTensor<2, Real> T2(values);

    T = 0;
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(0., T._coords[i], 1e-12);
    T2.zero();
    for (const auto i : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(0., T2._coords[i], 1e-12);
  }

  void testNorm()
  {
    std::vector<Real> values = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    TypeNTensor<2, Real> T(values);
    Real sum_sq = (8 * 9 * 17 / 6.);
    LIBMESH_ASSERT_FP_EQUAL(sum_sq, T.norm_sq(), 1e-12);
    LIBMESH_ASSERT_FP_EQUAL(std::sqrt(sum_sq), T.norm(), 1e-12);
  }

  void testSlice()
  {
    std::vector<Real> values = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    TypeNTensor<2, Real> T(values);
    const auto s0 = T.slice(0);
    const auto s1 = T.slice(1);
    const auto s2 = T.slice(2);
    for (const auto i : make_range(3))
      LIBMESH_ASSERT_FP_EQUAL(values[i], s0._coords[i], 1e-12);
    for (const auto i : make_range(3))
      LIBMESH_ASSERT_FP_EQUAL(values[3 + i], s1._coords[i], 1e-12);
    for (const auto i : make_range(3))
      LIBMESH_ASSERT_FP_EQUAL(values[6 + i], s2._coords[i], 1e-12);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TypeNTensorTest);
