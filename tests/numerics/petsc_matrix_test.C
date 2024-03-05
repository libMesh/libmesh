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

  CPPUNIT_TEST(testPetscBinaryRead);
  // CPPUNIT_TEST(testPetscHDF5Read);

  CPPUNIT_TEST_SUITE_END();

  void testPetscBinaryRead()
  {
    auto mat_to_read = std::make_unique<PetscMatrix<Number>>(*my_comm);

    // Petsc binary formats depend on sizeof(PetscInt)
#if LIBMESH_DOF_ID_BYTES == 4
    mat_to_read->read_petsc_binary("matrices/geom_1_extraction_op.petsc32");

    CPPUNIT_ASSERT_EQUAL(mat_to_read->m(), dof_id_type(27));
    CPPUNIT_ASSERT_EQUAL(mat_to_read->n(), dof_id_type(27));

    // Our read_matlab partitioning doesn't necessarily match PETSc's
    // MatLoad partitioning, so we'll partition a matrix to compare
    // against manually.  These particular files have bandwidth 8,
    // so we'll ask for sufficient n_nz and n_oz to handle that
    // regardless of partitioning.
    auto mat_ascii_format = std::make_unique<PetscMatrix<Number>>(*my_comm);
    mat_ascii_format->init(mat_to_read->m(), mat_to_read->n(),
                           mat_to_read->local_m(), mat_to_read->local_n(),
                           8, 7);

    mat_ascii_format->read_matlab("matrices/geom_1_extraction_op.m");

    mat_ascii_format->add(-1, *mat_to_read);
    CPPUNIT_ASSERT_LESS(TOLERANCE, mat_ascii_format->l1_norm());
#endif
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(PetscMatrixTest);

#endif
