// libmesh includes
#include <libmesh/lumped_mass_matrix.h>
#include <libmesh/parallel.h>
#include <libmesh/id_types.h>
#include <libmesh/int_range.h>
#include <libmesh/numeric_vector.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#include <memory>

#define UNUSED 0

using namespace libMesh;

class LumpedMassMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(LumpedMassMatrixTest);

  CPPUNIT_TEST(testNumerics);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
    _comm = TestCommWorld;
    _matrix = std::make_unique<LumpedMassMatrix<Number>>(*_comm);

    numeric_index_type root_block_size = 2;
    _local_size = root_block_size + static_cast<numeric_index_type>(_comm->rank());
    _global_size = 0;

    for (processor_id_type p = 0; p < _comm->size(); ++p)
    {
      numeric_index_type block_size = root_block_size + static_cast<numeric_index_type>(p);
      _global_size += block_size;
    }

    _matrix->init(_global_size, UNUSED, _local_size, UNUSED);
  }

  void tearDown() {}

  void testNumerics()
  {
    LOG_UNIT_TEST;

    numeric_index_type beginning_index = _matrix->row_start();
    numeric_index_type end_index = _matrix->row_stop();

    CPPUNIT_ASSERT_EQUAL(_local_size, numeric_index_type(_matrix->row_stop() - _matrix->row_start()));

    _matrix->zero();

    std::vector<Real> gold_values(_local_size, 0);

    // Test add
    for (const auto i : make_range(beginning_index, end_index))
      for (const auto j : make_range(beginning_index, end_index))
      {
        const Real sgn = j % 2 ? 1 : -1;
        _matrix->add(i, j, sgn * j);
        gold_values[i - beginning_index] += j;
      }

    _matrix->close();

    for (const auto i : make_range(beginning_index, end_index))
      for (const auto j : make_range(beginning_index, end_index))
      {
        if (i != j)
          LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real((*_matrix)(i, j)), _tolerance);
        else
          LIBMESH_ASSERT_FP_EQUAL(gold_values[i - beginning_index], libmesh_real((*_matrix)(i, i)), _tolerance);
      }

    // Test set
    for (const auto i : make_range(beginning_index, end_index))
    {
      const Real sgn = i % 2 ? 1 : -1;
      _matrix->set(i, i, sgn * gold_values[i - beginning_index]);
    }

    _matrix->close();

    for (const auto i : make_range(beginning_index, end_index))
      LIBMESH_ASSERT_FP_EQUAL(gold_values[i - beginning_index], libmesh_real((*_matrix)(i, i)), _tolerance);

    // Test clone
    auto copy = _matrix->clone();
    // Check that matrices have the same local/global sizes
    CPPUNIT_ASSERT_EQUAL(copy->m(), _matrix->m());
    CPPUNIT_ASSERT_EQUAL(copy->n(), _matrix->n());
    CPPUNIT_ASSERT_EQUAL(copy->local_m(), _matrix->local_m());
    CPPUNIT_ASSERT_EQUAL(copy->row_start(), _matrix->row_start());
    CPPUNIT_ASSERT_EQUAL(copy->row_stop(), _matrix->row_stop());

    for (const auto i : make_range(beginning_index, end_index))
      for (const auto j : make_range(beginning_index, end_index))
        LIBMESH_ASSERT_FP_EQUAL(libmesh_real((*_matrix)(i, j)), libmesh_real((*copy)(i, j)), _tolerance);

    // Test zero clone
    auto zero_copy = _matrix->zero_clone();
    // Check that matrices have the same local/global sizes
    CPPUNIT_ASSERT_EQUAL(zero_copy->m(), _matrix->m());
    CPPUNIT_ASSERT_EQUAL(zero_copy->n(), _matrix->n());
    CPPUNIT_ASSERT_EQUAL(zero_copy->local_m(), _matrix->local_m());
    CPPUNIT_ASSERT_EQUAL(zero_copy->row_start(), _matrix->row_start());
    CPPUNIT_ASSERT_EQUAL(zero_copy->row_stop(), _matrix->row_stop());

    for (const auto i : make_range(beginning_index, end_index))
      for (const auto j : make_range(beginning_index, end_index))
        LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real((*zero_copy)(i, j)), _tolerance);
  }

private:

  Parallel::Communicator * _comm;
  std::unique_ptr<LumpedMassMatrix<Number>> _matrix;
  numeric_index_type _local_size, _global_size;
  const Real _tolerance = TOLERANCE * TOLERANCE;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LumpedMassMatrixTest);
