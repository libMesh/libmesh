#include <libmesh/trilinos_epetra_vector.h>

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

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

class EpetraVectorTest : public NumericVectorTest<EpetraVector<Real> > {
public:
  CPPUNIT_TEST_SUITE( EpetraVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( EpetraVectorTest );

#endif // LIBMESH_TRILINOS_HAVE_EPETRA
