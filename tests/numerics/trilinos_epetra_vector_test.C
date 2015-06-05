#include <libmesh/trilinos_epetra_vector.h>

#ifdef LIBMESH_HAVE_TRILINOS

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

using namespace libMesh;

class EpetraVectorTest : public NumericVectorTest<EpetraVector<Real> > {
public:
  CPPUNIT_TEST_SUITE( EpetraVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( EpetraVectorTest );

#endif // #ifdef LIBMESH_HAVE_TRILINOS

