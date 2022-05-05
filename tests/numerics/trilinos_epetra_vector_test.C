#include <libmesh/trilinos_epetra_vector.h>

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

#include "numeric_vector_test.h"


using namespace libMesh;

class EpetraVectorTest : public NumericVectorTest<EpetraVector<Real>> {
public:
  EpetraVectorTest() :
    NumericVectorTest<EpetraVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "EpetraVectorTest";
  }

  CPPUNIT_TEST_SUITE( EpetraVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( EpetraVectorTest );

#endif // LIBMESH_TRILINOS_HAVE_EPETRA
