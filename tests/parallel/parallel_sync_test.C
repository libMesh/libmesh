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

  CPPUNIT_TEST( testPush );

//  CPPUNIT_TEST( testPushOversized );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testPushImpl(int M)
  {
    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();

    // Data to send/recieve with each processor rank.  For this test,
    // processor p will send to destination d the integer d, in a
    // vector with sqrt(c)+1 copies, iff c := |p-d| is a square number.
    std::map<processor_id_type, std::vector<unsigned int> > data, received_data;

    for (int d=0; d != M; ++d)
      {
        int diffsize = std::abs(d-rank);
        int diffsqrt = std::sqrt(diffsize);
        if (diffsqrt*diffsqrt == diffsize)
          for (int i=-1; i != diffsqrt; ++i)
            data[d].push_back(d);
      }

    auto collect_data =
      [&received_data]
      (processor_id_type pid,
       const std::vector<unsigned int> & data)
      {
        received_data[pid] = data;
      };

    // Do the push
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
    const int size = TestCommWorld->size();

    // Our sync functions are most typically used with a map of
    // processor ids that *only* includes ranks currently running.
    const int M = size;

    testPushImpl(M);
  }


  void testPushOversized()
  {
    const int size = TestCommWorld->size();

    // Our sync functions need to support sending to ranks that don't
    // exist!  If we're on N processors but working on a mesh
    // partitioned into M parts with M > N, then subpartition p
    // belongs to processor p%N.  Let's make M > N here.
    const int M = (size + 4) * 2;

    testPushImpl(M);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelSyncTest );
