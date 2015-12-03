#include <libmesh/laspack_vector.h>

#ifdef LIBMESH_HAVE_LASPACK

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

using namespace libMesh;

class LaspackVectorTest : public NumericVectorTest<LaspackVector<Real> > {
public:
  void setUp()
  {
    // Laspack doesn't support distributed parallel vectors, but we
    // can build a serial vector on each processor
    my_comm = new Parallel::Communicator();
  }

  void tearDown()
  {
    delete my_comm;
  }

  CPPUNIT_TEST_SUITE( LaspackVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( LaspackVectorTest );

#endif // #ifdef LIBMESH_HAVE_LASPACK

