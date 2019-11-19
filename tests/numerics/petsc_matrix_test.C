#include <libmesh/petsc_matrix.h>

#ifdef LIBMESH_HAVE_PETSC

#include <libmesh/parallel.h>
#include <libmesh/dense_matrix.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#include <vector>

using namespace libMesh;

class PetscMatrixTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(PetscMatrixTest);

  CPPUNIT_TEST(testGetAndSet);

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

    std::vector<Number> values;
    for (; count < _local_size; ++count, ++index)
    {
      _matrix->get_row(index, cols, values);
      for (numeric_index_type j = 0; j < _local_size; ++j)
      {
        LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(values[j]),
                                (count + 1) * (j + 1) * (_comm->rank() + 1),
                                _tolerance);
        CPPUNIT_ASSERT_EQUAL(cols[count], index);
      }
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
