// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/elem.h>
#include <libmesh/reference_elem.h>

// Unit test headers
#include "stream_redirector.h"

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

class WhichNodeAmITest : public CppUnit::TestCase
{

public:
  CPPUNIT_TEST_SUITE( WhichNodeAmITest );
  CPPUNIT_TEST( testPyramids );
  CPPUNIT_TEST( testPrisms );
  CPPUNIT_TEST( testTets );
  CPPUNIT_TEST( testHexes );
  CPPUNIT_TEST_SUITE_END();

public:

  void testPyramids()
  {
    // The last node on the right side (1) should be node 4 (apex node).
    const Elem & pyr5 = ReferenceElem::get(PYRAMID5);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4), pyr5.which_node_am_i(/*side=*/1, /*node=*/2));

    // Test the libmesh_asserts when they are enabled and exceptions
    // are available. If exceptions aren't available, libmesh_assert
    // simply aborts, so we can't unit test in that case.
#if !defined(NDEBUG) && defined(LIBMESH_ENABLE_EXCEPTIONS)
    try
      {
        // Avoid sending confusing error messages to the console.
        StreamRedirector stream_redirector;

        // Asking for the 4th node on a triangular face should throw.
        unsigned int n = pyr5.which_node_am_i(1, 3);

        // We shouldn't get here if the line above throws. If we do
        // get here, there's no way this assert will pass.
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(-1), n);
      }
    catch (...) {}
#endif

#ifdef NDEBUG
    // In optimized mode, we expect to get the "dummy" value 99.
    unsigned int n = pyr5.which_node_am_i(1, 3);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(99), n);
#endif

    // The last node on the right side (1) should be node 10.
    const Elem & pyr13 = ReferenceElem::get(PYRAMID13);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(10), pyr13.which_node_am_i(/*side=*/1, /*node=*/5));

    // The central node of the base should be node 13
    const Elem & pyr14 = ReferenceElem::get(PYRAMID14);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(13), pyr14.which_node_am_i(/*side=*/4, /*node=*/8));
  }



  void testPrisms()
  {
    // A PRISM6 has four nodes on some sides and three nodes on others
    const Elem & prism6 = ReferenceElem::get(PRISM6);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4), prism6.which_node_am_i(/*side=*/4, /*node=*/1));

    // Test the libmesh_asserts when they are enabled and exceptions
    // are available. If exceptions aren't available, libmesh_assert
    // simply aborts, so we can't unit test in that case.
#if !defined(NDEBUG) && defined(LIBMESH_ENABLE_EXCEPTIONS)
    try
      {
        // Avoid sending confusing error messages to the console.
        StreamRedirector stream_redirector;

        // Asks for the 3rd node on a Tri face which only has
        // indices 0, 1, and 2. Should throw an exception
        // (libmesh_assert throws an exception) when NDEBUG is not
        // defined.
        unsigned int n = prism6.which_node_am_i(0, 3);

        // We shouldn't get here if the line above throws. If we do
        // get here, there's no way this assert will pass.
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(-1), n);
      }
    catch (...) {}
#endif

#ifdef NDEBUG
    // In optimized mode, we expect to get the "dummy" value 99.
    unsigned int n = prism6.which_node_am_i(0, 3);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(99), n);
#endif

    // Test the Prism15.
    const Elem & prism15 = ReferenceElem::get(PRISM15);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), prism15.which_node_am_i(/*side=*/1, /*node=*/3));
  }



  void testTets()
  {
    const Elem & tet4 = ReferenceElem::get(TET4);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0), tet4.which_node_am_i(/*side=*/0, /*node=*/0));

    // Node 4 is a mid-edge node on side 1.
    const Elem & tet10 = ReferenceElem::get(TET10);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4), tet10.which_node_am_i(/*side=*/1, /*node=*/3));
  }



  void testHexes()
  {
    // Top left node on back side.
    const Elem & hex8 = ReferenceElem::get(HEX8);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(7), hex8.which_node_am_i(/*side=*/3, /*node=*/2));

    const Elem & hex20 = ReferenceElem::get(HEX20);
    const Elem & hex27 = ReferenceElem::get(HEX27);

    // The vertices (i.e. the first 4 nodes on each side should match for all of the Hex types.
    for (unsigned int side=0; side<hex8.n_sides(); ++side)
      for (unsigned int node=0; node<4; ++node)
        {
          // Make sure the Hex8 and Hex20 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex8.which_node_am_i(side, node),
                               hex20.which_node_am_i(side, node));

          // Make sure the Hex20 and Hex27 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex20.which_node_am_i(side, node),
                               hex27.which_node_am_i(side, node));
        }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( WhichNodeAmITest );
