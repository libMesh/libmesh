#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_LASPACK

#include <libmesh/laspack_matrix.h>

#include "sparse_matrix_test.h"

using namespace libMesh;

class LaspackMatrixTest : public SparseMatrixTest<LaspackMatrix<Number>>
{
public:
  LaspackMatrixTest() :
    SparseMatrixTest<LaspackMatrix<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "SparseMatrixTest";
    else
      this->libmesh_suite_name = "LaspackMatrixTest";
  }

  void setUp()
  {
    // LaspackMatrix is serial; we'll tell it to use MPI_COMM_SELF
    // so we just do these tests embarrassingly parallel
    my_comm = &comm_self;

    // LaspackMatrix doesn't support non-square matrices?
    nonsquare = 0;

    SparseMatrixTest<LaspackMatrix<Number>>::setUp();
  }

  CPPUNIT_TEST_SUITE(LaspackMatrixTest);

  SPARSEMATRIXTEST

  CPPUNIT_TEST_SUITE_END();

private:

  Parallel::Communicator comm_self;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LaspackMatrixTest);

#endif
