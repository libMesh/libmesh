#include <libmesh/eigen_sparse_vector.h>

#ifdef LIBMESH_HAVE_EIGEN

#include "numeric_vector_test.h"


using namespace libMesh;

class EigenSparseVectorTest : public NumericVectorTest<EigenSparseVector<libMesh::Number>> {
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

  EigenSparseVectorTest() :
    NumericVectorTest<EigenSparseVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "EigenSparseVectorTest";
  }

  CPPUNIT_TEST_SUITE( EigenSparseVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( EigenSparseVectorTest );

#endif // #ifdef LIBMESH_HAVE_EIGEN
