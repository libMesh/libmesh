#include <libmesh/vector_value.h>

#include "type_vector_test.h"

#include <type_traits>


using namespace libMesh;

#define VECTORVALUETEST                         \
  TYPEVECTORTEST                                \
  CPPUNIT_TEST( testScalarInit );               \

class RealVectorValueTest : public TypeVectorTestBase<VectorValue<Real>> {
public:
  CPPUNIT_TEST_SUITE( RealVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

class NumberVectorValueTest : public TypeVectorTestBase<VectorValue<Number>> {
public:
  CPPUNIT_TEST_SUITE( NumberVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

class ComplexVectorValueTest : public TypeVectorTestBase<VectorValue<Complex>> {
public:
  CPPUNIT_TEST_SUITE( NumberVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

class VectorCompareTypesTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( VectorCompareTypesTest );

  CPPUNIT_TEST( testCompareTypes );

  CPPUNIT_TEST_SUITE_END();

private:
  template <typename T, typename T2>
  typename CompareTypes<T, T2>::supertype
  average(const T & val1, const T2 & val2)
    {
      return (val1 + val2) / 2;
    }

public:
  void
  testCompareTypes()
    {
      VectorValue<float> fvec;
      VectorValue<double> dvec;

      auto ftype = fvec * 1;
      auto dtype = dvec * 1;

      {
        bool assertion = std::is_same<decltype(ftype), TypeVector<float>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        bool assertion = std::is_same<decltype(dtype), TypeVector<double>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(ftype, ftype);
        bool assertion = std::is_same<decltype(temp), TypeVector<float>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(ftype, dtype);
        bool assertion = std::is_same<decltype(temp), TypeVector<double>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(fvec, fvec);
        bool assertion = std::is_same<decltype(temp), VectorValue<float>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(fvec, dvec);
        bool assertion = std::is_same<decltype(temp), VectorValue<double>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(fvec, dtype);
        bool assertion = std::is_same<decltype(temp), VectorValue<double>>::value;
        CPPUNIT_ASSERT(assertion);
      }
      {
        auto temp = average(ftype, dvec);
        bool assertion = std::is_same<decltype(temp), VectorValue<double>>::value;
        CPPUNIT_ASSERT(assertion);
      }
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION( RealVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( NumberVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ComplexVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( VectorCompareTypesTest );
