// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel_sort.h>

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

class ParallelSortTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelSortTest );

  CPPUNIT_TEST( testSort );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testSort()
  {
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
