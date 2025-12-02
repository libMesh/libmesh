// libmesh includes
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

  CPPUNIT_TEST(testOperatorsScalar);
  CPPUNIT_TEST(testOperatorsTensor);
  CPPUNIT_TEST(testCastVector);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testSlice);

  CPPUNIT_TEST_SUITE_END();

private:
  void testOperatorsScalar()
  {
    LOG_UNIT_TEST;

    // Add a number

    // Subtract a number

    // Multiply by a number

    // Divide by a number

    // TensorValue<Real> tensor(1, 2, 0, 3, 4, 0);
    // VectorValue<Real> vector(5, 6, 0);
    // auto left_mult = vector * tensor;
    // auto right_mult = tensor * vector;
    // LIBMESH_ASSERT_FP_EQUAL(23, left_mult(0), 1e-12);
    // LIBMESH_ASSERT_FP_EQUAL(34, left_mult(1), 1e-12);
    // LIBMESH_ASSERT_FP_EQUAL(17, right_mult(0), 1e-12);
    // LIBMESH_ASSERT_FP_EQUAL(39, right_mult(1), 1e-12);
  }

  void testOperatorsTensor()
  {
    // Add a tensor

    // Return sum of two tensors

    // Subtract a tensor

    // Return subtraction of one tensor by another

    // Add a scaled tensor

  }

  void testCastVector()
  {

  }

  void testZero()
  {

  }

  void testSlice()
  {

  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TypeNTensorTest);
