#include <petsc_vector.h>

#ifdef HAVE_PETSC

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class PetscVectorTest : public NumericVectorTest<PetscVector<Real> > { 
public: 
  CPPUNIT_TEST_SUITE( PetscVectorTest );

  CPPUNIT_TEST( testLocalize );
  CPPUNIT_TEST( testLocalizeBase );
  CPPUNIT_TEST( testLocalizeToOne );
  CPPUNIT_TEST( testLocalizeToOneBase );

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( PetscVectorTest );

#endif // #ifdef HAVE_PETSC

