#include <libmesh/elem.h>
#include <libmesh/reference_elem.h>

// Unit test headers
#include "stream_redirector.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class WhichNodeAmITest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( WhichNodeAmITest );
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testPyramids );
  CPPUNIT_TEST( testPrisms );
  CPPUNIT_TEST( testTets );
  CPPUNIT_TEST( testHexes );
#endif
  CPPUNIT_TEST_SUITE_END();

public:

  void testPyramids()
  {
    LOG_UNIT_TEST;

    const Elem & pyr5 = ReferenceElem::get(PYRAMID5);
    // The last node on the right side (1) should be node 4 (apex node).
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4), pyr5.local_side_node(/*side=*/1, /*node=*/2));
    // (Edge 0, Node 0) -> Local node 0
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0), pyr5.local_edge_node(/*edge=*/0, /*node=*/0));
    // Second node of edges 4-7 is the apex node (4).
    for (unsigned int edge=4; edge<8; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4),
                           pyr5.local_edge_node(edge, /*node=*/1));

    // Test the libmesh_asserts when they are enabled and exceptions
    // are available. If exceptions aren't available, libmesh_assert
    // simply aborts, so we can't unit test in that case.
#if !defined(NDEBUG) && defined(LIBMESH_ENABLE_EXCEPTIONS)
    try
      {
        // Avoid sending confusing error messages to the console.
        StreamRedirector stream_redirector;

        // Asking for the 4th node on a triangular face should throw.
        unsigned int n = pyr5.local_side_node(1, 3);

        // We shouldn't get here if the line above throws. If we do
        // get here, there's no way this assert will pass.
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(-1), n);
      }
    catch (...) {}
#endif

#ifdef NDEBUG
    // In optimized mode, we expect to get the "dummy" value 99.
    unsigned int n = pyr5.local_side_node(1, 3);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(99), n);
#endif

    const Elem & pyr13 = ReferenceElem::get(PYRAMID13);
    // The last node on the right side (1) should be node 10.
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(10), pyr13.local_side_node(/*side=*/1, /*node=*/5));
    // Second node of edges 4-7 is the apex node (4).
    for (unsigned int edge=4; edge<8; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4),
                           pyr13.local_edge_node(edge, /*node=*/1));

    const Elem & pyr14 = ReferenceElem::get(PYRAMID14);
    // The central node of the base should be node 13
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(13), pyr14.local_side_node(/*side=*/4, /*node=*/8));
    // Second node of edges 4-7 is the apex node (4).
    for (unsigned int edge=4; edge<8; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4),
                           pyr14.local_edge_node(edge, /*node=*/1));

    const Elem & pyr18 = ReferenceElem::get(PYRAMID18);
    // The central node of the base should be node 13
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(13), pyr18.local_side_node(/*side=*/4, /*node=*/8));
    // The central node of the sides should be nodes 14-17
    for (unsigned int side=0; side<4; ++side)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(14+side), pyr18.local_side_node(side, /*node=*/6));
    // Second node of edges 4-7 is the apex node (4).
    for (unsigned int edge=4; edge<8; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4),
                           pyr18.local_edge_node(edge, /*node=*/1));
  }



  void testPrisms()
  {
    LOG_UNIT_TEST;

    // A PRISM6 has four nodes on some sides and three nodes on others
    for (ElemType elem_type : {PRISM6, PRISM15, PRISM18, PRISM20, PRISM21})
      {
        const Elem & prism = ReferenceElem::get(elem_type);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4),
                             prism.local_side_node(/*side=*/4, /*node=*/1));
        // Edges 3, 4, 5 are the "vertical" Prism edges.
        // Edges 0-2 and 6-8 are the "triangle" Prism edges.
        // We orient edges by volume-node order for some reason, which
        // makes some cases here weird.
        for (unsigned int i=0; i!=3; ++i)
          {
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i%2),
                                 prism.local_edge_node(/*edge=*/i, /*node=*/0));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1+(i>0)),
                                 prism.local_edge_node(/*edge=*/i, /*node=*/1));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i),
                                 prism.local_edge_node(/*edge=*/i+3, /*node=*/0));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+3),
                                 prism.local_edge_node(/*edge=*/i+3, /*node=*/1));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i%2+3),
                                 prism.local_edge_node(/*edge=*/i+6, /*node=*/0));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4+(i>0)),
                                 prism.local_edge_node(/*edge=*/i+6, /*node=*/1));
          }
      }

    // Test the libmesh_asserts when they are enabled and exceptions
    // are available. If exceptions aren't available, libmesh_assert
    // simply aborts, so we can't unit test in that case.
#if !defined(NDEBUG) && defined(LIBMESH_ENABLE_EXCEPTIONS)
    try
      {
        const Elem & prism6 = ReferenceElem::get(PRISM6);
        // Avoid sending confusing error messages to the console.
        StreamRedirector stream_redirector;

        // Asks for the 3rd node on a Tri face which only has
        // indices 0, 1, and 2. Should throw an exception
        // (libmesh_assert throws an exception) when NDEBUG is not
        // defined.
        unsigned int n = prism6.local_side_node(0, 3);

        // We shouldn't get here if the line above throws. If we do
        // get here, there's no way this assert will pass.
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(-1), n);
      }
    catch (...) {}
#endif

#ifdef NDEBUG
    // In optimized mode, we expect to get the "dummy" value 99.
    {
      const Elem & prism6 = ReferenceElem::get(PRISM6);
      unsigned int n = prism6.local_side_node(0, 3);
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(99), n);
    }
#endif

    // Test the Prism15 midedge nodes.
    for (ElemType elem_type : {PRISM15, PRISM18, PRISM20, PRISM21})
      {
        const Elem & prism = ReferenceElem::get(elem_type);
        for (unsigned int i=0; i!=3; ++i)
          {
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+6),
                                 prism.local_edge_node(/*edge=*/i, /*node=*/2));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+6),
                                 prism.local_side_node(/*side=*/i+1, /*node=*/4));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(8-i), // "flipped" side
                                 prism.local_side_node(/*side=*/0, /*node=*/i+3));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+9),
                                 prism.local_edge_node(/*edge=*/i+3, /*node=*/2));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+9),
                                 prism.local_side_node(/*side=*/i+1, /*node=*/7));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+9),
                                 prism.local_side_node(/*side=*/(i+2)%3+1, /*node=*/5));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+12),
                                 prism.local_edge_node(/*edge=*/i+6, /*node=*/2));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+12),
                                 prism.local_side_node(/*side=*/i+1, /*node=*/6));
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(i+12),
                                 prism.local_side_node(/*side=*/4, /*node=*/i+3));
          }
      }

    // Test the Prism18 face nodes.
    for (ElemType elem_type : {PRISM18, PRISM20, PRISM21})
      {
        const Elem & prism = ReferenceElem::get(elem_type);
        // Faces 1, 2, 3 are the "vertical" Prism sides
        for (unsigned int side=1; side<4; ++side)
          {
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(side+14),
                                 prism.local_side_node(/*side=*/side, /*node=*/8));
          }
      }

    // Test the Prism20 face nodes.
    for (ElemType elem_type : {PRISM20, PRISM21})
      {
        const Elem & prism = ReferenceElem::get(elem_type);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(18),
                             prism.local_side_node(/*side=*/0, /*node=*/6));
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(19),
                             prism.local_side_node(/*side=*/4, /*node=*/6));
      }

  }



  void testTets()
  {
    LOG_UNIT_TEST;

    const Elem & tet4 = ReferenceElem::get(TET4);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0), tet4.local_side_node(/*side=*/0, /*node=*/0));

    // Edges 3, 4, 5 all connect to the "apex" node 3
    for (unsigned int edge=3; edge<6; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3),
                           tet4.local_edge_node(edge, /*node=*/1));

    const Elem & tet10 = ReferenceElem::get(TET10);
    // Node 4 is a mid-edge node on side 1.
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(4), tet10.local_side_node(/*side=*/1, /*node=*/3));
    // Nodes 4, 5, 6 are mid-edge nodes of the first three edges
    for (unsigned int edge=0; edge<3; ++edge)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(edge + 4),
                           tet10.local_edge_node(edge, /*node=*/2));
  }



  void testHexes()
  {
    LOG_UNIT_TEST;

    // Top left node on back side.
    const Elem & hex8 = ReferenceElem::get(HEX8);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(7), hex8.local_side_node(/*side=*/3, /*node=*/2));

    const Elem & hex20 = ReferenceElem::get(HEX20);
    const Elem & hex27 = ReferenceElem::get(HEX27);

    // The vertices (i.e. the first 4 nodes on each side should match for all of the Hex types.
    for (unsigned int side=0; side<hex8.n_sides(); ++side)
      for (unsigned int node=0; node<4; ++node)
        {
          // Make sure the Hex8 and Hex20 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex8.local_side_node(side, node),
                               hex20.local_side_node(side, node));

          // Make sure the Hex20 and Hex27 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex20.local_side_node(side, node),
                               hex27.local_side_node(side, node));
        }

    for (unsigned int edge=0; edge<hex8.n_edges(); ++edge)
      for (unsigned int node=0; node<2; ++node)
        {
          // Make sure the Hex8 and Hex20 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex8.local_edge_node(edge, node),
                               hex20.local_edge_node(edge, node));

          // Make sure the Hex20 and Hex27 implementations agree.
          CPPUNIT_ASSERT_EQUAL(hex20.local_edge_node(edge, node),
                               hex27.local_edge_node(edge, node));
        }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( WhichNodeAmITest );
