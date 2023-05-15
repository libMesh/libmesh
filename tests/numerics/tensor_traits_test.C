#include "libmesh_cppunit.h"
#include <libmesh/tensor_tools.h>

using namespace libMesh;
using namespace TensorTools;

class TensorTraitsTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( TensorTraitsTest );

  CPPUNIT_TEST( test );

  CPPUNIT_TEST_SUITE_END();

public:
  void
  test()
    {
      LOG_UNIT_TEST;

      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(0), TensorTraits<Real>::rank);
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), TensorTraits<VectorValue<Real>>::rank);
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), TensorTraits<TypeVector<Real>>::rank);
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(2), TensorTraits<TensorValue<Real>>::rank);
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(2), TensorTraits<TypeTensor<Real>>::rank);
      typedef TypeNTensor<3, Real> TypeNTensorTestType;
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(3), TensorTraits<TypeNTensorTestType>::rank);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TensorTraitsTest );
