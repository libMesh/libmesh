#include <libmesh/parallel_sort.h>
#include <libmesh/parallel.h>

#include <algorithm>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class ParallelSortTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ParallelSortTest );

  CPPUNIT_TEST( testSort );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testSort()
  {
    LOG_UNIT_TEST;

    const int size = TestCommWorld->size(),
              rank = TestCommWorld->rank();
    const int n_vals = size - rank;
    std::vector<int> vals(n_vals);

    // Put all the numbers up to a triangular number into our data,
    // scrambled and backwards
    int val = rank+1, stride = size;
    for (int i=0; i != n_vals; ++i)
      {
        vals[n_vals-i-1] = val;
        val += stride;
        stride -= 1;
      }

    Parallel::Sort<int> sorter (*TestCommWorld, vals);

    sorter.sort();

    const std::vector<int> & my_bin = sorter.bin();

    // Our bins should be roughly the same size, but with that
    // nbins*50 stuff in Parallel::BinSorter it's hard to predict the
    // outcome exactly.  We'll just make sure they're sorted and
    // they've got everything.

    int total_size = cast_int<int>(my_bin.size());
    TestCommWorld->sum(total_size);

    CPPUNIT_ASSERT_EQUAL(total_size, size*(size+1)/2);

    CPPUNIT_ASSERT(std::is_sorted(my_bin.begin(), my_bin.end()));

    int rank_with_i = -1;
    for (int i=1; i <= total_size; ++i)
      {
        int count_i = std::count(my_bin.begin(), my_bin.end(), i);
        CPPUNIT_ASSERT(count_i < 2);

        if (count_i)
          {
            CPPUNIT_ASSERT(rank_with_i <= rank);
            rank_with_i = rank;
          }
        TestCommWorld->max(rank_with_i);

        TestCommWorld->sum(count_i);
        CPPUNIT_ASSERT_EQUAL(count_i, 1);
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelSortTest );
