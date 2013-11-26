// You can test libmesh's boost subset by seeing if this file compiles
// and runs against the boost installed in your $LIBMESH_DIR.
#include <boost/assert.hpp>
#include <boost/checked_delete.hpp>
#include <boost/current_function.hpp>
#include <boost/memory_order.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>
#include <boost/version.hpp>

// Workaround to avoid the following linking errors:
// Undefined symbols for architecture x86_64:
//   "boost::system::system_category()", referenced from:
//       __static_initialization_and_destruction_0(int, int) in cczaTREy.o
//   "boost::system::generic_category()", referenced from:
//       __static_initialization_and_destruction_0(int, int) in cczaTREy.o
// ld: symbol(s) not found for architecture x86_64
// collect2: error: ld returned 1 exit status
// which means that -lboost_system is required.
// See: http://stackoverflow.com/questions/17000542/boost-pool-can-i-wean-it-from-boost-system
#define BOOST_POOL_NO_MT       // disable multi-threading
#define BOOST_THREAD_MUTEX_HPP // define the #include-guard to disable the header

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/pool/object_pool.hpp>

#include <iostream>


// Helper function for testing the pool allocator
void test_pool();

int main (int, char**)
{
  boost::shared_ptr<int> p1( new int(1) );

  BOOST_ASSERT(true);

  std::cout << BOOST_CURRENT_FUNCTION << std::endl;

  test_pool();

  return 0;
}

void test_pool()
{
  boost::pool<> p(sizeof(int));
  for (int i = 0; i < 10000; ++i)
    {
      int * const t = static_cast<int* const>(p.malloc());
      *t = 42; // note: t is not freed
    }
} // on function exit, p is destroyed, and all malloc()'ed ints are implicitly freed
