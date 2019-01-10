// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel_sync.h>

#include "test_comm.h"

#include <algorithm>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

class ParallelSyncTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelSyncTest );

  // Our sync functions are most typically used with a map of
  // processor ids that *only* includes ranks currently running.
  CPPUNIT_TEST( testPush );
  CPPUNIT_TEST( testPull );
  CPPUNIT_TEST( testPushVecVec );
  CPPUNIT_TEST( testPullVecVec );

  // Our sync functions need to support sending to ranks that don't
  // exist!  If we're on N processors but working on a mesh
  // partitioned into M parts with M > N, then subpartition p belongs
  // to processor p%N.  Let's make M > N for these tests.
//  CPPUNIT_TEST( testPushOversized );
//  CPPUNIT_TEST( testPullOversized );
//  CPPUNIT_TEST( testPushVecVecOversized );
//  CPPUNIT_TEST( testPullVecVecOversized );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}


  // Data to send/recieve with each processor rank.  For this test,
  // processor p will send to destination d the integer d, in a vector
  // with sqrt(c)+1 copies, iff c := |p-d| is a square number.
  void fill_scalar_data
    (std::map<processor_id_type, std::vector<unsigned int>> & data,
     int M)
  {
    const int rank = TestCommWorld->rank();
    for (int d=0; d != M; ++d)
      {
        int diffsize = std::abs(d-rank);
        int diffsqrt = std::sqrt(diffsize);
        if (diffsqrt*diffsqrt == diffsize)
          for (int i=-1; i != diffsqrt; ++i)
            data[d].push_back(d);
      }
  }


  // Data to send/recieve with each processor rank.  For this test,
  // processor p will send to destination d the integer d, in two
  // subvectors with sqrt(c) and 1 copies, iff c := |p-d| is a square
  // number.
  void fill_vector_data
    (std::map<processor_id_type, std::vector<std::vector<unsigned int>>> & data,
     int M)
  {
    const int rank = TestCommWorld->rank();
    for (int d=0; d != M; ++d)
      {
        int diffsize = std::abs(d-rank);
        int diffsqrt = std::sqrt(diffsize);
        if (diffsqrt*diffsqrt == diffsize)
          {
            data[d].resize(2);
            for (int i=-1; i != diffsqrt; ++i)
              data[d][0].push_back(d);
            data[d][1].push_back(d);
          }
      }
  }

  template <typename MapType>
  std::function<void(processor_id_type, typename MapType::mapped_type)>
  data_collector (MapType & received_data)
  {
    return
      [&received_data]
      (processor_id_type pid,
       const typename MapType::mapped_type & data)
      {
        received_data[pid] = data;
      };
  }


  void testPushImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    std::map<processor_id_type, std::vector<unsigned int> > data, received_data;

    fill_scalar_data(data, M);

    Parallel::push_parallel_vector_data(*TestCommWorld, data, data_collector(received_data));

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          if (diffsqrt*diffsqrt != diffsize)
            {
              CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(0));
              continue;
            }

          CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(1));
          const std::vector<unsigned int> & datum = received_data[srcp];
          CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p), std::ptrdiff_t(diffsqrt+1));
          checked_sizes[srcp] += diffsqrt+1;
        }

    for (int srcp=0; srcp != size; ++srcp)
      CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], received_data[srcp].size());
  }


  void testPush()
  {
    testPushImpl(TestCommWorld->size());
  }


  void testPushOversized()
  {
    testPushImpl((TestCommWorld->size() + 4) * 2);
  }


  void testPullImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    std::map<processor_id_type, std::vector<unsigned int> > data, received_data;

    fill_scalar_data(data, M);

    auto compose_replies =
      []
      (processor_id_type pid,
       const std::vector<unsigned int> & query,
       std::vector<unsigned int> & response)
      {
        const std::size_t query_size = query.size();
        response.resize(query_size);
        for (unsigned int i=0; i != query_size; ++i)
          response[i] = query[i]*query[i];
      };


    auto collect_replies =
      [&received_data]
      (processor_id_type pid,
       const std::vector<unsigned int> & query,
       const std::vector<unsigned int> & response)
      {
        const std::size_t query_size = query.size();
        CPPUNIT_ASSERT_EQUAL(query_size, response.size());
        for (unsigned int i=0; i != query_size; ++i)
          {
            CPPUNIT_ASSERT_EQUAL(query[i]*query[i], response[i]);
          }
        received_data[pid] = response;
      };

    // Do the pull
    unsigned int * ex = nullptr;
    Parallel::pull_parallel_vector_data
      (*TestCommWorld, data, compose_replies, collect_replies, ex);

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          if (diffsqrt*diffsqrt != diffsize)
            {
              CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(0));
              continue;
            }

          CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(1));
          const std::vector<unsigned int> & datum = received_data[srcp];
          CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p*p), std::ptrdiff_t(diffsqrt+1));
          checked_sizes[srcp] += diffsqrt+1;
        }

    for (int srcp=0; srcp != size; ++srcp)
      CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], received_data[srcp].size());
  }


  void testPull()
  {
    testPullImpl(TestCommWorld->size());
  }


  void testPullOversized()
  {
    testPullImpl((TestCommWorld->size() + 4) * 2);
  }


  void testPushVecVecImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    std::map<processor_id_type, std::vector<std::vector<unsigned int>>> data, received_data;

    fill_vector_data(data, M);

    Parallel::push_parallel_vector_data(*TestCommWorld, data, data_collector(received_data));

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          if (diffsqrt*diffsqrt != diffsize)
            {
              CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(0));
              continue;
            }

          CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(1));
          CPPUNIT_ASSERT_EQUAL(received_data[srcp].size(), std::size_t(2));
          const std::vector<unsigned int> & datum = received_data[srcp][0];
          CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p), std::ptrdiff_t(diffsqrt+1));
          checked_sizes[srcp] += diffsqrt+1;
          CPPUNIT_ASSERT_EQUAL(received_data[srcp][1].size(), std::size_t(1));
          CPPUNIT_ASSERT_EQUAL(int(received_data[srcp][1][0]), p);
        }

    for (int srcp=0; srcp != size; ++srcp)
      CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], received_data[srcp][0].size());
  }


  void testPushVecVec()
  {
    testPushVecVecImpl(TestCommWorld->size());
  }


  void testPushVecVecOversized()
  {
    testPushVecVecImpl((TestCommWorld->size() + 4) * 2);
  }


  void testPullVecVecImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    std::map<processor_id_type, std::vector<std::vector<unsigned int>>> data, received_data;

    fill_vector_data(data, M);

    auto compose_replies =
      []
      (processor_id_type pid,
       const std::vector<std::vector<unsigned int>> & query,
       std::vector<std::vector<unsigned int>> & response)
      {
        const std::size_t query_size = query.size();
        response.resize(query_size);
        for (unsigned int i=0; i != query_size; ++i)
          {
            const std::size_t query_i_size = query[i].size();
            response[i].resize(query_i_size);
            for (unsigned int j=0; j != query_i_size; ++j)
            response[i][j] = query[i][j]*query[i][j];
          }
      };


    auto collect_replies =
      [&received_data]
      (processor_id_type pid,
       const std::vector<std::vector<unsigned int>> & query,
       const std::vector<std::vector<unsigned int>> & response)
      {
        const std::size_t query_size = query.size();
        CPPUNIT_ASSERT_EQUAL(query_size, response.size());
        for (unsigned int i=0; i != query_size; ++i)
          {
            const std::size_t query_i_size = query[i].size();
            CPPUNIT_ASSERT_EQUAL(query_i_size, response[i].size());
            for (unsigned int j=0; j != query_i_size; ++j)
              CPPUNIT_ASSERT_EQUAL(query[i][j]*query[i][j], response[i][j]);
          }
        received_data[pid] = response;
      };

    // Do the pull
    std::vector<unsigned int> * ex = nullptr;
    Parallel::pull_parallel_vector_data
      (*TestCommWorld, data, compose_replies, collect_replies, ex);

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          if (diffsqrt*diffsqrt != diffsize)
            {
              CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(0));
              continue;
            }

          CPPUNIT_ASSERT_EQUAL(received_data.count(srcp), std::size_t(1));
          CPPUNIT_ASSERT_EQUAL(received_data[srcp].size(), std::size_t(2));
          const std::vector<unsigned int> & datum = received_data[srcp][0];
          CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p*p), std::ptrdiff_t(diffsqrt+1));
          checked_sizes[srcp] += diffsqrt+1;
          CPPUNIT_ASSERT_EQUAL(received_data[srcp][1].size(), std::size_t(1));
          CPPUNIT_ASSERT_EQUAL(int(received_data[srcp][1][0]), p*p);
        }

    for (int srcp=0; srcp != size; ++srcp)
      CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], received_data[srcp][0].size());
  }


  void testPullVecVec()
  {
    testPullVecVecImpl(TestCommWorld->size());
  }


  void testPullVecVecOversized()
  {
    testPushVecVecImpl((TestCommWorld->size() + 4) * 2);
  }



};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelSyncTest );
