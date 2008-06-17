#include <distributed_vector.h>

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class DistributedVectorTest : public NumericVectorTest<DistributedVector<Real> > { 
public: 
  CPPUNIT_TEST_SUITE( DistributedVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( DistributedVectorTest );
