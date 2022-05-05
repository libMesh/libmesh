#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_EIGEN

// Unit test includes
#include "libmesh_cppunit.h"
#include "test_comm.h"

// libMesh includes
#include <libmesh/eigen_sparse_matrix.h>
#include <libmesh/dense_matrix.h>

// C++ includes
#include <memory>
#include <vector>

using namespace libMesh;

class EigenSparseMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(EigenSparseMatrixTest);

  CPPUNIT_TEST(testGetAndSet);
  CPPUNIT_TEST(testClone);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
    // EigenSparseMatrix is serial, but its constructor takes a comm
    // for consistency with other SparseMatrix types.
    _comm = TestCommWorld;
    _matrix = std::make_unique<EigenSparseMatrix<Number>>(*_comm);

    // All parameters are ignored except the number of global rows and colums and nnz.
    _matrix->init(/*m*/10,
                  /*n*/10,
                  /*m_l==m*/10,
                  /*n_l==n*/10,
                  /*nnz*/3);
  }

  void tearDown() {}

  void testGetAndSet()
  {
    LOG_UNIT_TEST;

    // EigenSparseMatrix is serial, so we simply test inserting the
    // same values on all procs.
    std::vector<numeric_index_type> rows = {0, 1, 2};
    std::vector<numeric_index_type> cols = {0, 1, 2};
    DenseMatrix<Number> local(3, 3);
    local.get_values() =
      {
        2., -1,  0,
        -1, 2., -1,
         0, -1,  2
      };

    _matrix->add_matrix(local, rows, cols);
    _matrix->close();
  }

  void testClone()
  {
    LOG_UNIT_TEST;

    {
      // Create copy, test that it can go out of scope
      auto copy = _matrix->clone();

      // Check that matrices have the same local/global sizes
      CPPUNIT_ASSERT_EQUAL(copy->m(), _matrix->m());
      CPPUNIT_ASSERT_EQUAL(copy->n(), _matrix->n());
      CPPUNIT_ASSERT_EQUAL(copy->local_m(), _matrix->local_m());
      CPPUNIT_ASSERT_EQUAL(copy->row_start(), _matrix->row_start());
      CPPUNIT_ASSERT_EQUAL(copy->row_stop(), _matrix->row_stop());

      // Check that copy has same values as original
      LIBMESH_ASSERT_FP_EQUAL(copy->l1_norm(), _matrix->l1_norm(), _tolerance);
    }

    {
      // Create zero copy
      auto zero_copy = _matrix->zero_clone();

      // Check that matrices have the same local/global sizes
      CPPUNIT_ASSERT_EQUAL(zero_copy->m(), _matrix->m());
      CPPUNIT_ASSERT_EQUAL(zero_copy->n(), _matrix->n());
      CPPUNIT_ASSERT_EQUAL(zero_copy->local_m(), _matrix->local_m());
      CPPUNIT_ASSERT_EQUAL(zero_copy->row_start(), _matrix->row_start());
      CPPUNIT_ASSERT_EQUAL(zero_copy->row_stop(), _matrix->row_stop());

      // Check that zero_copy has same values as original
      LIBMESH_ASSERT_FP_EQUAL(0.0, zero_copy->l1_norm(), _tolerance);
    }
  }

private:

  Parallel::Communicator * _comm;
  std::unique_ptr<EigenSparseMatrix<Number>> _matrix;
  const Real _tolerance = TOLERANCE * TOLERANCE;

};

CPPUNIT_TEST_SUITE_REGISTRATION(EigenSparseMatrixTest);

#endif
