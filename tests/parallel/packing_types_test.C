#include <libmesh/parallel.h>
#include <libmesh/parallel_eigen.h>
#include <libmesh/int_range.h>

#include "timpi/parallel_sync.h"

#include <vector>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

#ifdef LIBMESH_HAVE_EIGEN

class PackingTypesTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( PackingTypesTest );

  CPPUNIT_TEST( testDynamicEigenMatrix );
  CPPUNIT_TEST( testDynamicEigenVector );
  CPPUNIT_TEST( testNonFixedScalar );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



  void testDynamicEigenMatrix()
  {
    LOG_UNIT_TEST;

    typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> DynamicEigenMatrix;
    typedef std::vector<DynamicEigenMatrix> DynamicMatrixVector;
    std::map<processor_id_type, DynamicMatrixVector> send_data;

    // prepare unique data
    for (const auto i : index_range(*TestCommWorld))
    {
      const auto j = TestCommWorld->rank();
      send_data[i] = DynamicMatrixVector{ DynamicEigenMatrix(i+1, i+2) };
      for (const auto row: make_range(i+1))
        for (const auto col: make_range(i+2))
          send_data[i][0](row, col) = row * 100000.0 + col * 1000.0 + 10.0 * i + j;
    }

    // verification
    auto verify_data =
    [](processor_id_type j, const DynamicMatrixVector & recv_data)
    {
      const auto i = TestCommWorld->rank();
      const std::size_t rows = recv_data[0].rows();
      const std::size_t cols = recv_data[0].cols();

      CPPUNIT_ASSERT_EQUAL( rows, static_cast<std::size_t>(i+1) );
      CPPUNIT_ASSERT_EQUAL( cols, static_cast<std::size_t>(i+2) );

      for (const auto row: make_range(rows))
        for (const auto col: make_range(cols))
          CPPUNIT_ASSERT_EQUAL( recv_data[0](row, col), static_cast<Real>(row * 100000.0 + col * 1000.0 + j + 10.0 * i) );
    };

    // communicate
    Parallel::push_parallel_vector_data(*TestCommWorld, send_data, verify_data);
  }

  void testDynamicEigenVector()
  {
    LOG_UNIT_TEST;

    typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> DynamicEigenVector;
    typedef std::vector<DynamicEigenVector> DynamicVectorVector;
    std::map<processor_id_type, DynamicVectorVector> send_data;

    // prepare unique data
    for (const auto i : index_range(*TestCommWorld))
    {
      const auto j = TestCommWorld->rank();
      send_data[i] = DynamicVectorVector{ DynamicEigenVector(i+1, 1) };
      for (const auto row: make_range(i+1))
          send_data[i][0](row) = row * 1000.0 + 10.0 * i + j;
    }

    // verification
    auto verify_data =
    [](processor_id_type j, const DynamicVectorVector & recv_data)
    {
      const auto i = TestCommWorld->rank();
      const std::size_t rows = recv_data[0].rows();
      const std::size_t cols = recv_data[0].cols();

      CPPUNIT_ASSERT_EQUAL( rows, static_cast<std::size_t>(i+1) );
      CPPUNIT_ASSERT_EQUAL( cols, 1lu );

      for (const auto row: make_range(rows))
        CPPUNIT_ASSERT_EQUAL( recv_data[0](row), static_cast<Real>(row * 1000.0 + j + 10.0 * i) );
    };

    // communicate
    Parallel::push_parallel_vector_data(*TestCommWorld, send_data, verify_data);
  }

  void testNonFixedScalar()
  {
    LOG_UNIT_TEST;

    typedef Eigen::Matrix<std::vector<Real>, 3, 3> VectorEigenMatrix;
    typedef std::vector<VectorEigenMatrix> VectorMatrixVector;
    std::map<processor_id_type, VectorMatrixVector> send_data;

    // prepare unique data
    for (const auto i : index_range(*TestCommWorld))
    {
      const auto j = TestCommWorld->rank();
      send_data[i] = VectorMatrixVector(1);
      for (const auto row: make_range(3))
        for (const auto col: make_range(3))
          send_data[i][0](row, col).assign(row + col + 1, 10.0 * i + j);
    }

    // verification
    auto verify_data =
    [](processor_id_type j, const VectorMatrixVector & recv_data)
    {
      const auto i = TestCommWorld->rank();
      const std::size_t rows = recv_data[0].rows();
      const std::size_t cols = recv_data[0].cols();

      CPPUNIT_ASSERT_EQUAL( rows, static_cast<std::size_t>(3) );
      CPPUNIT_ASSERT_EQUAL( cols, static_cast<std::size_t>(3) );

      for (const auto row: make_range(rows))
        for (const auto col: make_range(cols))
        {
          CPPUNIT_ASSERT_EQUAL( recv_data[0](row, col).size(), static_cast<std::size_t>(row + col + 1) );
          for (const auto val : recv_data[0](row, col))
            CPPUNIT_ASSERT_EQUAL( val, static_cast<Real>(j + 10.0 * i) );
        }
    };

    // communicate
    Parallel::push_parallel_vector_data(*TestCommWorld, send_data, verify_data);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( PackingTypesTest );

#endif // LIBMESH_HAVE_EIGEN
