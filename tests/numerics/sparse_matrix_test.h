#ifndef SPARSE_MATRIX_TEST_H
#define SPARSE_MATRIX_TEST_H

// test includes
#include "test_comm.h"

// libMesh includes
#include <libmesh/dense_matrix.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/fuzzy_equals.h>

#include "libmesh_cppunit.h"

// C++ includes
#include <vector>
#ifdef LIBMESH_HAVE_CXX11_THREAD
#include <thread>
#include <algorithm>
#endif


#define SPARSEMATRIXTEST                     \
  CPPUNIT_TEST(testGetAndSet);               \
  CPPUNIT_TEST(testReadMatlab1);             \
  CPPUNIT_TEST(testReadMatlab2);             \
  CPPUNIT_TEST(testReadMatlab4);             \
  CPPUNIT_TEST(testWriteAndRead);            \
  CPPUNIT_TEST(testClone);



template <class DerivedClass>
class SparseMatrixTest : public CppUnit::TestCase
{
public:

  void setUp()
  {
    // By default we'll use the whole communicator in parallel;
    // Serial-only NumericVector subclasses will need to override
    // this and set something else first.
    if (!my_comm)
      my_comm = TestCommWorld;

    matrix = std::make_unique<DerivedClass>(*my_comm);

    // Use the same even partitioning that we'll auto-deduce in matrix
    // read, to make it easier to test parallel reads

    global_m = global_n = my_comm->size() * 5.75 - 1;
    if (nonsquare)
      global_n *= 1.3;

    libMesh::dof_id_type
      row_start =     my_comm->rank() * global_m / my_comm->size(),
      row_stop  = (my_comm->rank()+1) * global_m / my_comm->size();
    libMesh::dof_id_type
      col_start =     my_comm->rank() * global_n / my_comm->size(),
      col_stop  = (my_comm->rank()+1) * global_n / my_comm->size();

    local_m = row_stop - row_start;
    local_n = col_stop - col_start;

    // Let's just play around with locally dense blocks for now
    matrix->init(global_m,
                 global_n,
                 local_m,
                 local_n,
                 /*nnz=*/local_n,
                 /*noz=*/5);
  }


  void tearDown() {}


  void setValues()
  {
    std::vector<libMesh::numeric_index_type> rows(local_m);
    std::iota(rows.begin(), rows.end(), matrix->row_start());
    std::vector<libMesh::numeric_index_type> cols(local_n);
    std::iota(cols.begin(), cols.end(), matrix->col_start());

    libMesh::DenseMatrix<libMesh::Number> local(local_m, local_n);

    for (auto i : libMesh::make_range(local_m))
      for (auto j : libMesh::make_range(local_n))
        local(i, j) = (i + 1) * (j + 1) * (my_comm->rank() + 1);

    matrix->zero();

    matrix->add_matrix(local, rows, cols);
    matrix->close();
  }


  void testValues()
  {
    auto functor = [this]()
    {
      std::vector<libMesh::numeric_index_type> cols_to_get;
      std::vector<libMesh::Number> values;
      const libMesh::numeric_index_type col_start =
        matrix->col_start();
      for (libMesh::numeric_index_type i :
           libMesh::make_range(matrix->row_start(),
                               matrix->row_stop()))
        {
          matrix->get_row(i, cols_to_get, values);
          for (libMesh::numeric_index_type col_j :
               libMesh::index_range(cols_to_get))
            {
              CPPUNIT_ASSERT_EQUAL(cols_to_get[col_j], col_start + col_j);
              LIBMESH_ASSERT_FP_EQUAL
                ((i - matrix->row_start() + 1) * (col_j + 1) * (my_comm->rank() + 1),
                 libMesh::libmesh_real(values[col_j]), _tolerance);
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


  void testGetAndSet()
  {
    LOG_UNIT_TEST;

    setValues();

    testValues();
  }

  void testReadMatlab(const std::string & filename)
  {
    // Laspack doesn't handle non-square matrices)
    if (matrix->solver_package() == libMesh::LASPACK_SOLVERS)
      return;

    matrix->clear();

    auto matrix2 = std::make_unique<DerivedClass>(*my_comm);

    matrix->read(filename);

#ifndef LIBMESH_HAVE_GZSTREAM
    return;
#endif

    matrix2->read(std::string(filename)+".gz");

    // We need some more SparseMatrix operators, but not today
    CPPUNIT_ASSERT(matrix->l1_norm() == matrix2->l1_norm());
    CPPUNIT_ASSERT(matrix->linfty_norm() == matrix2->linfty_norm());
  }


  void testReadMatlab1()
  {
    LOG_UNIT_TEST;
    testReadMatlab("matrices/geom_1_extraction_op.m");
  }


  void testReadMatlab2()
  {
    LOG_UNIT_TEST;
    testReadMatlab("matrices/geom_2_extraction_op.m");
  }


  void testReadMatlab4()
  {
    LOG_UNIT_TEST;
    testReadMatlab("matrices/geom_4_extraction_op.m");
  }


  void testWriteAndRead()
  {
    LOG_UNIT_TEST;

    setValues();

    // If we're working with serial matrices then just print one of
    // them so they don't step on the others' toes.
    if (matrix->n_processors() > 1 ||
        TestCommWorld->rank() == 0)
      matrix->print_matlab(libmesh_suite_name+"_matrix.m");

    matrix->clear();

    // Let's make sure we don't have any race conditions; we have
    // multiple SparseMatrix subclasses that might be trying to read
    // and write the same file.

    TestCommWorld->barrier();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    // We're not supporting complex reads quite yet
    return;
#endif

    matrix->read(libmesh_suite_name+"_matrix.m");

    TestCommWorld->barrier();

    testValues();
  }

  void testClone()
  {
    LOG_UNIT_TEST;

    setValues();

    // Matrix must be closed before it can be cloned.
    matrix->close();

    {
      // Create copy, test that it can go out of scope
      auto copy = matrix->clone();

      // Check that matrices have the same local/global sizes
      CPPUNIT_ASSERT_EQUAL(copy->m(), matrix->m());
      CPPUNIT_ASSERT_EQUAL(copy->n(), matrix->n());
      CPPUNIT_ASSERT_EQUAL(copy->local_m(), matrix->local_m());
      CPPUNIT_ASSERT_EQUAL(copy->row_start(), matrix->row_start());
      CPPUNIT_ASSERT_EQUAL(copy->row_stop(), matrix->row_stop());

      // Check that copy has same values as original
      LIBMESH_ASSERT_FP_EQUAL(copy->l1_norm(), matrix->l1_norm(), _tolerance);
      CPPUNIT_ASSERT(relative_fuzzy_equals(*matrix, *copy));
      copy->scale(2);
      CPPUNIT_ASSERT(!relative_fuzzy_equals(*matrix, *copy));
    }

    {
      // Create zero copy
      auto zero_copy = matrix->zero_clone();

      // Check that matrices have the same local/global sizes
      CPPUNIT_ASSERT_EQUAL(zero_copy->m(), matrix->m());
      CPPUNIT_ASSERT_EQUAL(zero_copy->n(), matrix->n());
      CPPUNIT_ASSERT_EQUAL(zero_copy->local_m(), matrix->local_m());
      CPPUNIT_ASSERT_EQUAL(zero_copy->row_start(), matrix->row_start());
      CPPUNIT_ASSERT_EQUAL(zero_copy->row_stop(), matrix->row_stop());

      // Check that zero_copy has same values as original
      LIBMESH_ASSERT_FP_EQUAL(0.0, zero_copy->l1_norm(), _tolerance);
    }
  }

protected:

  std::string libmesh_suite_name;

  libMesh::Parallel::Communicator * my_comm = nullptr;

  std::unique_ptr<DerivedClass> matrix;

  libMesh::numeric_index_type nonsquare = 1,
      local_m, local_n,
      global_m, global_n;

  const libMesh::Real _tolerance = libMesh::TOLERANCE * libMesh::TOLERANCE;
};

#endif // SPARSE_MATRIX_TEST_H
