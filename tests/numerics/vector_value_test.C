// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/vector_value.h>

#include "type_vector_test.h"

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

CPPUNIT_TEST_SUITE_REGISTRATION( RealVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( NumberVectorValueTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ComplexVectorValueTest );
