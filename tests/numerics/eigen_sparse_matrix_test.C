#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_EIGEN

#include <libmesh/eigen_sparse_matrix.h>

#include "sparse_matrix_test.h"

using namespace libMesh;

class EigenSparseMatrixTest : public SparseMatrixTest<EigenSparseMatrix<Number>>
{
public:
  EigenSparseMatrixTest() :
    SparseMatrixTest<EigenSparseMatrix<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "SparseMatrixTest";
    else
      this->libmesh_suite_name = "EigenSparseMatrixTest";
  }

  void setUp()
  {
    // EigenSparseMatrix is serial; we'll tell it to use MPI_COMM_SELF
    // so we just do these tests embarrassingly parallel
    my_comm = &comm_self;

    // EigenSparseMatrix doesn't support non-square matrices?
    nonsquare = 0;

    SparseMatrixTest<EigenSparseMatrix<Number>>::setUp();
  }

  CPPUNIT_TEST_SUITE(EigenSparseMatrixTest);

  SPARSEMATRIXTEST

  CPPUNIT_TEST_SUITE_END();

private:

  Parallel::Communicator comm_self;
};

CPPUNIT_TEST_SUITE_REGISTRATION(EigenSparseMatrixTest);

#endif
