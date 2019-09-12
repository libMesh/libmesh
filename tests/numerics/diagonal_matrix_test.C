// libmesh includes
#include <libmesh/diagonal_matrix.h>
#include <libmesh/parallel.h>
#include <libmesh/auto_ptr.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/id_types.h>
#include <libmesh/numeric_vector.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#include <memory>
#include <numeric>

#define UNUSED 0

using namespace libMesh;

namespace
{
template <typename T>
T
termial(T n)
{
  if (n == 1)
    return 1;
  return n + termial(n - 1);
}
}

class DiagonalMatrixTest : public CppUnit::TestCase
{
public:
  void setUp()
  {
    _comm = TestCommWorld;
    _matrix = libmesh_make_unique<DiagonalMatrix<Real>>(*_comm);

    numeric_index_type root_block_size = 2;
    _local_size = root_block_size + static_cast<numeric_index_type>(_comm->rank());
    _global_size = 0;
    _i.push_back(0);

    for (processor_id_type p = 0; p < _comm->size(); ++p)
    {
      numeric_index_type block_size = root_block_size + static_cast<numeric_index_type>(p);
      _global_size += block_size;
      _i.push_back(block_size);
    }

    _matrix->init(_global_size, UNUSED, _local_size, UNUSED);
  }

  void tearDown() {}

  CPPUNIT_TEST_SUITE(DiagonalMatrixTest);

  CPPUNIT_TEST(testSizes);
  CPPUNIT_TEST(testNumerics);

  CPPUNIT_TEST_SUITE_END();

private:
  void testSizes()
  {
    CPPUNIT_ASSERT_EQUAL(_global_size, _matrix->m());
    CPPUNIT_ASSERT_EQUAL(_global_size, _matrix->n());
    CPPUNIT_ASSERT_EQUAL(_i[_comm->rank()], _matrix->row_start());
    CPPUNIT_ASSERT_EQUAL(_i[_comm->rank() + 1], _matrix->row_stop());
  }

  void testNumerics()
  {
    numeric_index_type beginning_index = _matrix->row_start();
    numeric_index_type end_index = _matrix->row_stop();

    CPPUNIT_ASSERT_EQUAL(_local_size, _matrix->row_stop() - _matrix->row_start());

    _matrix->zero();

    // In general we don't want to mix these calls without an intervening close(), but since these
    // calls are not actually doing anything for this matrix type it's ok here
    _matrix->add(beginning_index, beginning_index + 1, 1);
    _matrix->set(beginning_index + 1, beginning_index, 1);

    // Test that off-diagonal elements are always zero
    LIBMESH_ASSERT_FP_EQUAL((*_matrix)(beginning_index, beginning_index + 1), 0, _tolerance);
    LIBMESH_ASSERT_FP_EQUAL((*_matrix)(beginning_index + 1, beginning_index), 0, _tolerance);

    // Test add
    {
      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        _matrix->add(i, i, i);

      // Also add to rank 0
      _matrix->add(0, 0, 1);

      CPPUNIT_ASSERT(!_matrix->closed());

      _matrix->close();
      CPPUNIT_ASSERT(_matrix->closed());

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
      {
        if (i)
          LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), i, _tolerance);
        else
          LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), _comm->size(), _tolerance);
      }
    }

    // Test set
    {
      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        if (i != _global_size - 1)
          _matrix->set(i, i, i);

      if (_comm->rank() == 0)
        _matrix->set(_global_size - 1, _global_size - 1, 0);

      CPPUNIT_ASSERT(!_matrix->closed());

      _matrix->close();
      CPPUNIT_ASSERT(_matrix->closed());

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
      {
        if (i != _global_size - 1)
          LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), i, _tolerance);
        else
          LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), 0, _tolerance);
      }
    }

    std::vector<numeric_index_type> rows(_local_size);
    std::iota(rows.begin(), rows.end(), beginning_index);

    // Test dense matrix add, get diagonal, SparseMatrix add, get transpose
    {
      DenseMatrix<Real> dense(_local_size, _local_size);

      for (numeric_index_type local_index = 0, global_index = beginning_index;
           local_index < _local_size;
           ++local_index, ++global_index)
        dense(local_index, local_index) = global_index;

      _matrix->zero();

      // rows-columns overload for DenseMatrix add
      _matrix->add_matrix(dense, rows, rows);

      // Close not really necessary here because we didn't do anything off-process but still good
      // practice
      _matrix->close();

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), i, _tolerance);

      // Build
      auto diagonal = NumericVector<Real>::build(*_comm);
      // Allocate storage
      diagonal->init(_matrix->diagonal());
      // Fill entries
      _matrix->get_diagonal(*diagonal);

      // Build
      DiagonalMatrix<Real> copy(*_comm);
      // Allocate storage
      copy.init(*_matrix);
      // Fill entries from diagonal
      copy = std::move(*diagonal);

      // Single dof-indices overload for DenseMatrix add
      _matrix->add_matrix(dense, rows);
      _matrix->close();

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), 2 * i, _tolerance);

      // SparseMatrix add
      _matrix->add(-1, copy);

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, i), i, _tolerance);

      _matrix->get_transpose(copy);

      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        LIBMESH_ASSERT_FP_EQUAL(copy(i, i), i, _tolerance);
    }

    // Test norms
    {
      LIBMESH_ASSERT_FP_EQUAL(termial(_global_size - 1), _matrix->l1_norm(), _tolerance);
      LIBMESH_ASSERT_FP_EQUAL(_global_size - 1, _matrix->linfty_norm(), _tolerance);
    }

    // Test zero_rows
    {
      _matrix->zero_rows(rows, 1);
      for (numeric_index_type i = beginning_index; i < end_index; ++i)
        for (numeric_index_type j = beginning_index; j < end_index; ++j)
        {
          if (i == j)
            LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, j), 1, _tolerance);
          else
            LIBMESH_ASSERT_FP_EQUAL((*_matrix)(i, j), 0, _tolerance);
        }
    }
  }

  Parallel::Communicator * _comm;
  std::unique_ptr<DiagonalMatrix<Real>> _matrix;
  numeric_index_type _local_size, _global_size;
  std::vector<numeric_index_type> _i;
  const Real _tolerance = TOLERANCE * TOLERANCE;
};

CPPUNIT_TEST_SUITE_REGISTRATION(DiagonalMatrixTest);
