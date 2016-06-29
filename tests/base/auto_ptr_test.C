// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <iomanip>
#include <libmesh/auto_ptr.h>

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

class AutoPtrTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( AutoPtrTest );

  // Test that we can call if (!foo) for AutoPtr
  CPPUNIT_TEST( testComparison );

  // Test that the libmesh_make_unique macro works.
  CPPUNIT_TEST( testMakeUnique );

  CPPUNIT_TEST_SUITE_END();

  // Note: this "extra" public declaration needs to be here, otherwise
  // there is an unmatched private section coming from one of the
  // macros above.
public:
  void setUp ()
  {}

  void tearDown ()
  {}

  void testComparison ()
  {
    // Test that if(foo) and if (!foo) compile and do the right thing
    // for AutoPtr.  This tests that the safe_bool thing that AutoPtr
    // uses is actually working.
    {
      AutoPtr<int> foo(new int(42));

      bool test1_passed = false;
      if (foo)
        test1_passed = true;
      CPPUNIT_ASSERT(test1_passed);

      bool test2_passed = true;
      if (!foo)
        test2_passed = false;
      CPPUNIT_ASSERT(test2_passed);
    }

    // Test the converse for when foo holds a NULL pointer.
    {
      AutoPtr<int> foo(NULL);

      bool test3_passed = true;
      if (foo)
        test3_passed = false;
      CPPUNIT_ASSERT(test3_passed);

      bool test4_passed = false;
      if (!foo)
        test4_passed = true;
      CPPUNIT_ASSERT(test4_passed);
    }

    // Make sure that foo == bar is not allowed.  This should not even
    // compile -- I don't think that's possible to test with CPPUnit
    // so it's commented out for now, but if you uncomment this test,
    // you should get a compiler error.  This is probably a candidate
    // for a configure-time test instead.
    // AutoPtr<int> bar(new int(21));
    // if (foo == bar)
    //   std::cerr << "This should not compile." << std::endl;
  }

  void testMakeUnique ()
  {
    // If the compiler supports one of the ways of calling
    // make_unique, test that our macro works for calling it works.
#if defined(LIBMESH_HAVE_CXX14_MAKE_UNIQUE) || defined(LIBMESH_HAVE_CXX11_MAKE_UNIQUE_WORKAROUND)
    std::unique_ptr<int> foo = libmesh_make_unique<int>(42);
#endif
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( AutoPtrTest );
