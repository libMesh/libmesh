#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_PETSC

#include <libmesh/petsc_matrix.h>

#include "sparse_matrix_test.h"

using namespace libMesh;

class PetscMatrixTest : public SparseMatrixTest<PetscMatrix<Number>>
{
public:
  PetscMatrixTest() :
    SparseMatrixTest<PetscMatrix<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "SparseMatrixTest";
    else
      this->libmesh_suite_name = "PetscMatrixTest";
  }

  CPPUNIT_TEST_SUITE(PetscMatrixTest);

  SPARSEMATRIXTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION(PetscMatrixTest);

#endif
