// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/fe_lagrange_shape_1D.h"

// unit test includes
#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

/**
 * This class is for unit testing lagrange shape function values and derivatives
 */
class LagrangeShapeTest : public CppUnit::TestCase
{
public:

  LIBMESH_CPPUNIT_TEST_SUITE( LagrangeShapeTest );
  CPPUNIT_TEST( test1DLagrange );
  CPPUNIT_TEST_SUITE_END();

public:

  void test1DLagrange ()
  {
    LOG_UNIT_TEST;

    Real tol = TOLERANCE*TOLERANCE;

    for (auto xi : {0., 1./3., -.5, -1./7.})
      for(auto o : make_range(1, 4))
        for(auto i : make_range(o + 1))
        {
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape(Order(o), i, xi),
                                  fe_lagrange_1D_any_shape(Order(o), i, xi), tol);
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape_deriv(Order(o), i, 0, xi),
                                  fe_lagrange_1D_any_shape_deriv(Order(o), i, 0, xi), tol);
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape_second_deriv(Order(o), i, 0, xi),
                                  fe_lagrange_1D_any_shape_second_deriv(Order(o), i, 0, xi), tol);
        }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( LagrangeShapeTest );
