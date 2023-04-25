#include <libmesh/eigen_sparse_vector.h>

#ifdef LIBMESH_HAVE_EIGEN

#include "numeric_vector_test.h"


using namespace libMesh;

class EigenSparseVectorTest : public NumericVectorTest<EigenSparseVector<Number>> {
private:
  // This class manages the memory for the communicator it uses,
  // providing a dumb pointer to the managed resource for the base
  // class.
  std::unique_ptr<Parallel::Communicator> _managed_comm;

public:
  void setUp()
  {
    // Eigen doesn't support distributed parallel vectors, but we can
    // build a serial vector on each processor
    _managed_comm = std::make_unique<Parallel::Communicator>();

    // Base class communicator points to our managed communicator
    my_comm = _managed_comm.get();

    this->NumericVectorTest<EigenSparseVector<Number>>::setUp();
  }

  void tearDown() {}

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
