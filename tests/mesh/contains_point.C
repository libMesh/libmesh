// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh.h>
#include <libmesh/elem.h>

#include "test_comm.h"

using namespace libMesh;

class ContainsPointTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the contains_point
   * method that is used extensively in the point locator. This test focusses on
   * the specializes contains_point implementaion in TRI3.
   */
public:
  CPPUNIT_TEST_SUITE( ContainsPointTest );

  CPPUNIT_TEST( testContainsPointTri3 );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  // TRI3 test
  void testContainsPointTri3() {
    // tri corner nodes
    Node a(3,1,2, 0);
    Node b(1,2,3, 1);
    Node c(2,3,1, 2);

    // helper vectors to span the tri and point out of plane
    Point va(a-c);
    Point vb(b-c);
    Point oop(va.cross(vb));

    UniquePtr<Elem> elem = Elem::build(TRI3);

    elem->set_node(0) = &a;
    elem->set_node(1) = &b;
    elem->set_node(2) = &c;

    // midpoint
    Point mid = 1.0/3.0 * (a + b + c);
    CPPUNIT_ASSERT (elem->contains_point(mid));

    // out of plane from the mid point
    CPPUNIT_ASSERT (!elem->contains_point(mid + 0.1 * oop));

    // check all corners
    CPPUNIT_ASSERT (elem->contains_point(a));
    CPPUNIT_ASSERT (elem->contains_point(b));
    CPPUNIT_ASSERT (elem->contains_point(c));

    // check outside
    CPPUNIT_ASSERT (!elem->contains_point(a + va * TOLERANCE * 10));
    CPPUNIT_ASSERT (!elem->contains_point(b + vb * TOLERANCE * 10));
    CPPUNIT_ASSERT (!elem->contains_point(c - (va + vb) * TOLERANCE * 10));
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( ContainsPointTest );
