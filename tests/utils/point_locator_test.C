// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/node.h>

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



class PointLocatorTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( PointLocatorTest );

  CPPUNIT_TEST( testLocatorOnEdge3 );
  CPPUNIT_TEST( testLocatorOnQuad9 );
  CPPUNIT_TEST( testLocatorOnTri6 );
  CPPUNIT_TEST( testLocatorOnHex27 );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testLocator(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    const unsigned n_elem_per_side = 5;
    const UniquePtr<Elem> test_elem = Elem::build(elem_type);
    const Real ymax = test_elem->dim() > 1;
    const Real zmax = test_elem->dim() > 2;
    const unsigned int ny = ymax * n_elem_per_side;
    const unsigned int nz = zmax * n_elem_per_side;

    MeshTools::Generation::build_cube (mesh,
                                       n_elem_per_side,
                                       ny,
                                       nz,
                                       0., 1.,
                                       0., ymax,
                                       0., zmax,
                                       elem_type);

    UniquePtr<PointLocatorBase> locator = mesh.sub_point_locator();

    if (!mesh.is_serial())
      locator->enable_out_of_mesh_mode();

    for (unsigned int i=0; i != n_elem_per_side+1; ++i)
      {
        for (unsigned int j=0; j != ny+1; ++j)
          {
            for (unsigned int k=0; k != nz+1; ++k)
              {
                const libMesh::Real h = libMesh::Real(1)/n_elem_per_side;
                Point p(i*h, j*h, k*h);

                const Elem *elem = locator->operator()(p);

                bool found_elem = elem;
                if (!mesh.is_serial())
                  mesh.comm().max(found_elem);

                CPPUNIT_ASSERT(found_elem);
                if (elem)
                  {
                    CPPUNIT_ASSERT(elem->contains_point(p));
                  }

                const Node *node = locator->locate_node(p);

                bool found_node = node;
                if (!mesh.is_serial())
                  mesh.comm().max(found_node);

                CPPUNIT_ASSERT(found_node);

                if (node)
                  {
                    CPPUNIT_ASSERT_DOUBLES_EQUAL((*node)(0), i*h,
                                                 TOLERANCE*TOLERANCE);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL((*node)(1), j*h,
                                                 TOLERANCE*TOLERANCE);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL((*node)(2), k*h,
                                                 TOLERANCE*TOLERANCE);
                  }
              }
          }
      }
  }



  void testLocatorOnEdge3() { testLocator(EDGE3); }
  void testLocatorOnQuad9() { testLocator(QUAD9); }
  void testLocatorOnTri6()  { testLocator(TRI6); }
  void testLocatorOnHex27() { testLocator(HEX27); }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PointLocatorTest );
