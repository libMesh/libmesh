// libMesh includes
#include "libmesh/utility.h"

// cppunit includes
#include "libmesh_cppunit.h"

// C++ includes
#include <set>
#include <memory>

using namespace libMesh;

class TransparentComparatorTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE ( TransparentComparatorTest );
  CPPUNIT_TEST( testSet );
  CPPUNIT_TEST_SUITE_END();

  void testSet()
  {
    LOG_UNIT_TEST;

    std::set<std::unique_ptr<int>, Utility::CompareUnderlying> s;
    s.insert(std::make_unique<int>(1));
    s.insert(std::make_unique<int>(2));
    s.insert(std::make_unique<int>(3));
    s.insert(std::make_unique<int>(3)); // Not a duplicate. This set does pointer comparisons, not value comparisons.

    // Search for entry by pointer-to-int. This tests the desired
    // behavior of the std::set comparison object.
    int * first = s.begin()->get();
    auto result = s.find(first);
    CPPUNIT_ASSERT(result != s.end());
    CPPUNIT_ASSERT(result->get() == first);

    // Test set::count(). It should be using the same underlying
    // machinery as std::set(find) but it's good to verify this.
    CPPUNIT_ASSERT(s.count(first) == 1);

    // Test that attempting to insert the same underlying pointer again fails. We simulate this by
    // creating and moving from a standalone object.
    auto four = std::make_unique<int>(4);

    // The first insert() adds a non-nullptr underlying pointer to the
    // set and leaves "four" with and underlying nullptr.
    s.insert(std::move(four));
    CPPUNIT_ASSERT(s.size() == 5);

    // The second insert() adds a nullptr underlying pointer to the
    // set and leaves "four" with and underlying nullptr.
    s.insert(std::move(four));
    CPPUNIT_ASSERT(s.size() == 6);

    // The third insert() does not change the set because it tries to
    // add a second nullptr underlying pointer to the set.
    s.insert(std::move(four));
    CPPUNIT_ASSERT(s.size() == 6);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION ( TransparentComparatorTest );
