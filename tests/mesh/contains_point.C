// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh.h>
#include <libmesh/elem.h>

#include "test_comm.h"

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

class ContainsPointTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the contains_point
   * method that is used extensively in the point locator. This test focuses on
   * the specializes contains_point implementation in TRI3.
   */
public:
  CPPUNIT_TEST_SUITE( ContainsPointTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testContainsPointTri3 );
  CPPUNIT_TEST( testContainsPointTet4 );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  // TRI3 test
  void testContainsPointTri3()
  {
    Point a(3,1,2), b(1,2,3), c(2,3,1);
    containsPointTri3Helper(a, b, c, a);

    // Ben's 1st "small triangle" test case.
    containsPointTri3Helper(a/5000., b/5000., c/5000., c/5000.);

    // Ben's 2nd "small triangle" test case.
    containsPointTri3Helper(Point(0.000808805, 0.0047284,  0.),
                            Point(0.0011453,   0.00472831, 0.),
                            Point(0.000982249, 0.00508037, 0.),
                            Point(0.001,       0.005,      0.));
  }



  // TET4 test
  void testContainsPointTet4()
  {
    // Construct unit Tet.
    {
      Node zero  (0., 0., 0., 0);
      Node one   (1., 0., 0., 1);
      Node two   (0., 1., 0., 2);
      Node three (0., 0., 1., 3);

      std::unique_ptr<Elem> elem = Elem::build(TET4);
      elem->set_node(0) = &zero;
      elem->set_node(1) = &one;
      elem->set_node(2) = &two;
      elem->set_node(3) = &three;

      // The centroid must be inside the element.
      CPPUNIT_ASSERT (elem->contains_point(elem->centroid()));

      // The vertices must be contained in the element.
      CPPUNIT_ASSERT (elem->contains_point(zero));
      CPPUNIT_ASSERT (elem->contains_point(one));
      CPPUNIT_ASSERT (elem->contains_point(two));
      CPPUNIT_ASSERT (elem->contains_point(three));

      // Make sure that outside points are not contained.
      CPPUNIT_ASSERT (!elem->contains_point(Point(.34, .34, .34)));
      CPPUNIT_ASSERT (!elem->contains_point(Point(.33, .33, -.1)));
      CPPUNIT_ASSERT (!elem->contains_point(Point(0., -.1, .5)));
    }


    // Construct a nearly degenerate (sliver) tet.  A unit tet with
    // nodes "one" and "two" moved to within a distance of epsilon
    // from the origin.
    {
      Real epsilon = 1.e-4;

      Node zero  (0., 0., 0., 0);
      Node one   (epsilon, 0., 0., 1);
      Node two   (0., epsilon, 0., 2);
      Node three (0., 0., 1., 3);

      std::unique_ptr<Elem> elem = Elem::build(TET4);
      elem->set_node(0) = &zero;
      elem->set_node(1) = &one;
      elem->set_node(2) = &two;
      elem->set_node(3) = &three;

      // The centroid must be inside the element.
      CPPUNIT_ASSERT (elem->contains_point(elem->centroid()));

      // The vertices must be contained in the element.
      CPPUNIT_ASSERT (elem->contains_point(zero));
      CPPUNIT_ASSERT (elem->contains_point(one));
      CPPUNIT_ASSERT (elem->contains_point(two));
      CPPUNIT_ASSERT (elem->contains_point(three));

      // Verify that a point which is on a mid-edge is contained in the element.
      CPPUNIT_ASSERT (elem->contains_point(Point(epsilon/2, 0, 0.5)));

      // Make sure that "just barely" outside points are outside.
      CPPUNIT_ASSERT (!elem->contains_point(Point(epsilon, epsilon, epsilon/2)));
      CPPUNIT_ASSERT (!elem->contains_point(Point(epsilon/10, epsilon/10, 1.0)));
      CPPUNIT_ASSERT (!elem->contains_point(Point(epsilon/2, -epsilon/10, 0.5)));
    }
  }

protected:
  // The first 3 arguments are the points of the triangle, the last argument is a
  // custom point that should be inside the Tri.
  void containsPointTri3Helper(Point a_in, Point b_in, Point c_in, Point p)
  {
    // vertices
    Node a(a_in, 0);
    Node b(b_in, 1);
    Node c(c_in, 2);

    // Build the test Elem
    std::unique_ptr<Elem> elem = Elem::build(TRI3);

    elem->set_node(0) = &a;
    elem->set_node(1) = &b;
    elem->set_node(2) = &c;

    // helper vectors to span the tri and point out of plane
    Point va(a-c);
    Point vb(b-c);
    Point oop(va.cross(vb));
    Point oop_unit = oop.unit();

    // The triangle must contain its own centroid
    CPPUNIT_ASSERT (elem->contains_point(elem->centroid()));

    // triangle should contain all its vertices
    CPPUNIT_ASSERT (elem->contains_point(a));
    CPPUNIT_ASSERT (elem->contains_point(b));
    CPPUNIT_ASSERT (elem->contains_point(c));

    // out of plane from the centroid, should *not* be in the element
    CPPUNIT_ASSERT (!elem->contains_point(elem->centroid() + std::sqrt(TOLERANCE) * oop_unit));

    // These points should be outside the triangle
    CPPUNIT_ASSERT (!elem->contains_point(a + va * TOLERANCE * 10));
    CPPUNIT_ASSERT (!elem->contains_point(b + vb * TOLERANCE * 10));
    CPPUNIT_ASSERT (!elem->contains_point(c - (va + vb) * TOLERANCE * 10));

    // Test the custom point
    CPPUNIT_ASSERT (elem->contains_point(p));
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( ContainsPointTest );
