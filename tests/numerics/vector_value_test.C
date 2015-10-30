// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/vector_value.h>

#include "type_vector_test.h"

using namespace libMesh;

#define VECTORVALUETEST                         \
  TYPEVECTORTEST                                \
  CPPUNIT_TEST( testScalarInit );               \

class RealVectorValueTest : public TypeVectorTestBase<VectorValue<Real> > {
public:
  CPPUNIT_TEST_SUITE( RealVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

class NumberVectorValueTest : public TypeVectorTestBase<VectorValue<Number> > {
public:
  CPPUNIT_TEST_SUITE( NumberVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

class ComplexVectorValueTest : public TypeVectorTestBase<VectorValue<Complex> > {
public:
  CPPUNIT_TEST_SUITE( NumberVectorValueTest );

  VECTORVALUETEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( RealVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( NumberVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ComplexVectorValueTest );
