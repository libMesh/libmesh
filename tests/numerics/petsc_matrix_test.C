#include <libmesh/petsc_matrix.h>

#ifdef LIBMESH_HAVE_PETSC

#include <libmesh/parallel.h>
#include <libmesh/dense_matrix.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#include <vector>
#ifdef LIBMESH_HAVE_CXX11_THREAD
#include <thread>
#include <algorithm>
#endif

using namespace libMesh;

class PetscMatrixTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(PetscMatrixTest);

  CPPUNIT_TEST(testGetAndSet);
  CPPUNIT_TEST(testClone);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
    _comm = TestCommWorld;
    _matrix = libmesh_make_unique<PetscMatrix<Number>>(*_comm);

    numeric_index_type root_block_size = 2;
    _local_size = root_block_size + static_cast<numeric_index_type>(_comm->rank());
    _global_size = 0;
    _i.push_back(0);

    for (processor_id_type p = 0; p < _comm->size(); ++p)
    {
      numeric_index_type block_size = root_block_size + static_cast<numeric_index_type>(p);
      _global_size += block_size;
      _i.push_back(_global_size);
    }

    _matrix->init(_global_size,
                  _global_size,
                  _local_size,
                  _local_size,
                  /*nnz=*/_local_size,
                  /*noz=*/0);
  }

  void tearDown() {}

  void testGetAndSet()
  {
    std::vector<numeric_index_type> rows(_local_size);
    std::vector<numeric_index_type> cols(_local_size);
    DenseMatrix<Number> local(_local_size, _local_size);

    numeric_index_type index = _i[_comm->rank()], count = 0;
    for (; count < _local_size; ++count, ++index)
    {
      rows[count] = index;
      cols[count] = index;
      for (numeric_index_type j = 0; j < _local_size; ++j)
        local(count, j) = (count + 1) * (j + 1) * (_comm->rank() + 1);
    }

    _matrix->add_matrix(local, rows, cols);
    _matrix->close();

    index = _i[_comm->rank()], count = 0;

    auto functor = [this]()
    {
      std::vector<numeric_index_type> cols_to_get;
      std::vector<Number> values;
      numeric_index_type local_index = _i[_comm->rank()], local_count = 0;
      for (; local_count < _local_size; ++local_count, ++local_index)
      {
        _matrix->get_row(local_index, cols_to_get, values);
        for (numeric_index_type j = 0; j < _local_size; ++j)
        {
          LIBMESH_ASSERT_FP_EQUAL((local_count + 1) * (j + 1) * (_comm->rank() + 1),
                                  libMesh::libmesh_real(values[j]),
                                  _tolerance);
          CPPUNIT_ASSERT_EQUAL(cols_to_get[local_count], local_index);
        }
      }
    };

#ifdef LIBMESH_HAVE_CXX11_THREAD
    auto num_threads = std::min(unsigned(2),
                                std::max(
                                  std::thread::hardware_concurrency(),
                                  unsigned(1)));
    std::vector<std::thread> threads(num_threads);
    for (unsigned int thread = 0; thread < num_threads; ++thread)
      threads[thread] = std::thread(functor);
    std::for_each(threads.begin(), threads.end(),
                  [](std::thread & x){x.join();});
#else
    functor();
#endif
  }

  void testClone()
  {
    // Matrix must be closed before it can be cloned.
    _matrix->close();

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
  std::unique_ptr<PetscMatrix<Number>> _matrix;
  numeric_index_type _local_size, _global_size;
  std::vector<numeric_index_type> _i;
  const Real _tolerance = TOLERANCE * TOLERANCE;

};

CPPUNIT_TEST_SUITE_REGISTRATION(PetscMatrixTest);

#endif
