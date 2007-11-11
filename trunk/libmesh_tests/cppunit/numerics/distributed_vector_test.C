#include <distributed_vector.h>

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class DistributedVectorTest : public NumericVectorTest<DistributedVector<Real> > { 
public: 
  CPPUNIT_TEST_SUITE( DistributedVectorTest );

  CPPUNIT_TEST( testLocalize );
  CPPUNIT_TEST( testLocalizeBase );
  CPPUNIT_TEST( testLocalizeToOne );
  CPPUNIT_TEST( testLocalizeToOneBase );

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( DistributedVectorTest );
