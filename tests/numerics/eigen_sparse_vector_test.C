#include <libmesh/eigen_sparse_vector.h>

#ifdef LIBMESH_HAVE_EIGEN

#include "numeric_vector_test.h"

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

using namespace libMesh;

class EigenSparseVectorTest : public NumericVectorTest<EigenSparseVector<libMesh::Number> > {
public:
  void setUp()
  {
    // Eigen doesn't support distributed parallel vectors, but we can
    // build a serial vector on each processor
    my_comm = new Parallel::Communicator();
  }

  void tearDown()
  {
    delete my_comm;
  }

  CPPUNIT_TEST_SUITE( EigenSparseVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( EigenSparseVectorTest );

#endif // #ifdef LIBMESH_HAVE_EIGEN
