#include <laspack_vector.h>

#ifdef HAVE_LASPACK

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class LaspackVectorTest : public NumericVectorTest<LaspackVector<Real> > { 
public: 
  CPPUNIT_TEST_SUITE( LaspackVectorTest );

  NUMERICVECTORTEST
  
  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( LaspackVectorTest );

#endif // #ifdef HAVE_LASPACK

