// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel_sync.h>
#include <libmesh/int_range.h>
#include <libmesh/simple_range.h>

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
  CPPUNIT_TEST( testPushMultimap );
  CPPUNIT_TEST( testPushMultimapVecVec );

  // Our sync functions need to support sending to ranks that don't
  // exist!  If we're on N processors but working on a mesh
  // partitioned into M parts with M > N, then subpartition p belongs
  // to processor p%N.  Let's make M > N for these tests.
  CPPUNIT_TEST( testPushOversized );
  CPPUNIT_TEST( testPullOversized );
  CPPUNIT_TEST( testPushVecVecOversized );
  CPPUNIT_TEST( testPullVecVecOversized );
  CPPUNIT_TEST( testPushMultimapOversized );
  CPPUNIT_TEST( testPushMultimapVecVecOversized );

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


  // Multimap data to send/recieve with each processor rank.  For this
  // test, processor p will send to destination d the integer d, in a
  // vector with sqrt(c)+1 copies followed by a vector with 1 copy,
  // iff c := |p-d| is a square number.
  void fill_scalar_data
    (std::multimap<processor_id_type, std::vector<unsigned int>> & data,
     int M)
  {
    const int rank = TestCommWorld->rank();
    for (int d=0; d != M; ++d)
      {
        int diffsize = std::abs(d-rank);
        int diffsqrt = std::sqrt(diffsize);
        if (diffsqrt*diffsqrt == diffsize)
          {
            std::vector<unsigned int> v;
            for (int i=-1; i != diffsqrt; ++i)
              v.push_back(d);
            data.emplace(d, v);
            v.resize(1, d);
            data.emplace(d, v);
          }
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



  // Multimap data to send/recieve with each processor rank.  For this
  // test, processor p will send to destination d the integer d, in
  // two subvectors with sqrt(c) and 1 copies, followed by a vector
  // with 1 copy, iff c := |p-d| is a square number.
  void fill_vector_data
    (std::multimap<processor_id_type, std::vector<std::vector<unsigned int>>> & data,
     int M)
  {
    const int rank = TestCommWorld->rank();
    for (int d=0; d != M; ++d)
      {
        int diffsize = std::abs(d-rank);
        int diffsqrt = std::sqrt(diffsize);
        if (diffsqrt*diffsqrt == diffsize)
          {
            std::vector<std::vector<unsigned int>> vv(2);
            for (int i=-1; i != diffsqrt; ++i)
              vv[0].push_back(d);
            vv[1].push_back(d);
            data.emplace(d, vv);
            vv.resize(1);
            vv[0].resize(1);
            data.emplace(d, vv);
          }
      }
  }


  void testPushImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    std::map<processor_id_type, std::vector<unsigned int> > data, received_data;

    fill_scalar_data(data, M);

    auto collect_data =
      [&received_data]
      (processor_id_type pid,
       const typename std::vector<unsigned int> & data)
      {
        auto & vec = received_data[pid];
        vec.insert(vec.end(), data.begin(), data.end());
      };

    Parallel::push_parallel_vector_data(*TestCommWorld, data, collect_data);

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
              if (received_data.count(srcp))
                {
                  const std::vector<unsigned int> & datum = received_data[srcp];
                  CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p), std::ptrdiff_t(0));
                }
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

    // Test the received results, for each query we sent.
    for (int p=0; p != M; ++p)
      {
        CPPUNIT_ASSERT_EQUAL(data[p].size(), received_data[p].size());
        for (auto i : index_range(data[p]))
          CPPUNIT_ASSERT_EQUAL(data[p][i]*data[p][i], received_data[p][i]);
      }
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

    std::map<processor_id_type, std::vector<std::vector<unsigned int>>> data;
    std::map<processor_id_type, std::vector<unsigned int>> received_data;

    fill_vector_data(data, M);

    auto collect_data =
      [&received_data]
      (processor_id_type pid,
       const typename std::vector<std::vector<unsigned int>> & data)
      {
        auto & vec = received_data[pid];
        vec.insert(vec.end(), data[0].begin(), data[0].end());
        CPPUNIT_ASSERT_EQUAL(data.size(), std::size_t(2));
        CPPUNIT_ASSERT_EQUAL(data[1].size(), std::size_t(1));
        CPPUNIT_ASSERT_EQUAL(data[0][0], data[1][0]);
      };

    Parallel::push_parallel_vector_data(*TestCommWorld, data, collect_data);

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
              if (received_data.count(srcp))
                {
                  const std::vector<unsigned int> & datum = received_data[srcp];
                  CPPUNIT_ASSERT_EQUAL(std::count(datum.begin(), datum.end(), p), std::ptrdiff_t(0));
                }
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
    std::map<processor_id_type, std::vector<std::vector<unsigned int>>> data;
    std::map<processor_id_type, std::vector<std::vector<unsigned int>>> received_data;

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
        auto & vec = received_data[pid];
        vec.emplace_back(response[0].begin(), response[0].end());
        CPPUNIT_ASSERT_EQUAL(response[1].size(), std::size_t(1));
        CPPUNIT_ASSERT_EQUAL(response[1][0], response[0][0]);
        vec.emplace_back(response[1].begin(), response[1].end());
      };

    // Do the pull
    std::vector<unsigned int> * ex = nullptr;
    Parallel::pull_parallel_vector_data
      (*TestCommWorld, data, compose_replies, collect_replies, ex);

    // Test the received results, for each query we sent.
    for (int p=0; p != M; ++p)
      {
        CPPUNIT_ASSERT_EQUAL(data[p].size(), received_data[p].size());
        for (auto i : index_range(data[p]))
          for (auto j : index_range(data[p][i]))
            CPPUNIT_ASSERT_EQUAL(data[p][i][j]*data[p][i][j], received_data[p][i][j]);
      }
  }


  void testPullVecVec()
  {
    testPullVecVecImpl(TestCommWorld->size());
  }


  void testPullVecVecOversized()
  {
    testPushVecVecImpl((TestCommWorld->size() + 4) * 2);
  }


  void testPushMultimapImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    // This is going to make sense because of C++11's guarantees
    // regarding preservation of insert ordering in multimaps,
    // combined with MPI's guarantees about non-overtaking
    std::multimap<processor_id_type, std::vector<unsigned int> > data, received_data;

    fill_scalar_data(data, M);

    auto collect_data =
      [&received_data]
      (processor_id_type pid,
       const typename std::vector<unsigned int> & data)
      {
        received_data.emplace(pid, data);
      };

    Parallel::push_parallel_vector_data(*TestCommWorld, data, collect_data);

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          auto rng = received_data.equal_range(srcp);
          if (diffsqrt*diffsqrt != diffsize)
            {
              for (auto & pv_it : as_range(rng))
                {
                  CPPUNIT_ASSERT_EQUAL(std::count(pv_it.second.begin(), pv_it.second.end(), p), std::ptrdiff_t(0));
                }
              continue;
            }

          CPPUNIT_ASSERT(rng.first != rng.second);
          for (auto pv_it = rng.first; pv_it != rng.second; ++pv_it)
            {
              std::ptrdiff_t cnt = std::count(pv_it->second.begin(), pv_it->second.end(), p);
              if (cnt)
                {
                  CPPUNIT_ASSERT_EQUAL(cnt, std::ptrdiff_t(diffsqrt+1));
                  auto pv_it2 = pv_it; ++pv_it2;
                  CPPUNIT_ASSERT(pv_it2 != rng.second);
                  std::ptrdiff_t cnt2 = std::count(pv_it2->second.begin(), pv_it2->second.end(), p);
                  CPPUNIT_ASSERT_EQUAL(cnt2, std::ptrdiff_t(1));
                  checked_sizes[srcp] += cnt + cnt2;
                  break;
                }
            }
        }

    for (int srcp=0; srcp != size; ++srcp)
      {
        std::size_t total_size = 0;
        for (auto & pv_it : as_range(received_data.equal_range(srcp)))
          total_size += pv_it.second.size();
        CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], total_size);
      }
  }


  void testPushMultimap()
  {
    testPushMultimapImpl(TestCommWorld->size());
  }


  void testPushMultimapOversized()
  {
    testPushMultimapImpl((TestCommWorld->size() + 4) * 2);
  }


  void testPushMultimapVecVecImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    // This is going to make sense because of C++11's guarantees
    // regarding preservation of insert ordering in multimaps,
    // combined with MPI's guarantees about non-overtaking
    std::multimap<processor_id_type, std::vector<std::vector<unsigned int>>> data, received_data;

    fill_vector_data(data, M);

    auto collect_data =
      [&received_data]
      (processor_id_type pid,
       const typename std::vector<std::vector<unsigned int>> & data)
      {
        received_data.emplace(pid, data);
      };

    Parallel::push_parallel_vector_data(*TestCommWorld, data, collect_data);

    // Test the received results, for each processor id p we're in
    // charge of.
    std::vector<std::size_t> checked_sizes(size, 0);
    for (int p=rank; p != M; p += size)
      for (int srcp=0; srcp != size; ++srcp)
        {
          int diffsize = std::abs(srcp-p);
          int diffsqrt = std::sqrt(diffsize);
          auto rng = received_data.equal_range(srcp);
          if (diffsqrt*diffsqrt != diffsize)
            {
              for (auto & pvv : as_range(rng))
                {
                  for (auto & v : pvv.second)
                  CPPUNIT_ASSERT_EQUAL(std::count(v.begin(), v.end(), p), std::ptrdiff_t(0));
                }
              continue;
            }

          CPPUNIT_ASSERT(rng.first != rng.second);
          for (auto pvv_it = rng.first; pvv_it != rng.second; ++pvv_it)
            {
              if(pvv_it->second.size() != std::size_t(2))
                libmesh_error();
              CPPUNIT_ASSERT_EQUAL(pvv_it->second.size(), std::size_t(2));
              std::ptrdiff_t cnt = std::count(pvv_it->second[0].begin(), pvv_it->second[0].end(), p);
              if (cnt)
                {
                  CPPUNIT_ASSERT_EQUAL(cnt, std::ptrdiff_t(diffsqrt+1));
                  std::ptrdiff_t cnt2 = std::count(pvv_it->second[1].begin(), pvv_it->second[1].end(), p);
                  CPPUNIT_ASSERT_EQUAL(cnt2, std::ptrdiff_t(1));
                  auto pvv_it2 = pvv_it; ++pvv_it2;
                  CPPUNIT_ASSERT(pvv_it2 != rng.second);
                  CPPUNIT_ASSERT_EQUAL(pvv_it2->second.size(), std::size_t(1));
                  std::ptrdiff_t cnt3 = std::count(pvv_it2->second[0].begin(), pvv_it2->second[0].end(), p);
                  CPPUNIT_ASSERT_EQUAL(cnt3, std::ptrdiff_t(1));
                  checked_sizes[srcp] += cnt + cnt2 + cnt3;
                  break;
                }
              ++pvv_it;
              libmesh_assert(pvv_it != rng.second);
            }
        }

    for (int srcp=0; srcp != size; ++srcp)
      {
        std::size_t total_size = 0;
        for (auto & pvv : as_range(received_data.equal_range(srcp)))
          for (auto & v : pvv.second)
            total_size += v.size();
        CPPUNIT_ASSERT_EQUAL(checked_sizes[srcp], total_size);
      }
  }


  void testPushMultimapVecVec()
  {
    testPushMultimapVecVecImpl(TestCommWorld->size());
  }


  void testPushMultimapVecVecOversized()
  {
    testPushMultimapVecVecImpl((TestCommWorld->size() + 4) * 2);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelSyncTest );
