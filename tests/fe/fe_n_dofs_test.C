#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/face_c0polygon.h>
#include <libmesh/fe.h>

#include "libmesh_cppunit.h"

using namespace libMesh;

class FENDofsTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to pin down the dof counts that
   * fe_lagrange_n_dofs.h and fe_monomial_n_dofs.h state, including the
   * cases those headers cannot answer: an order a family does not have,
   * and a polytope, whose count follows its node count rather than its
   * element type.
   *
   * FE<Dim,FAMILY> has a comma in it, which does not survive as a macro
   * argument, hence the typedefs.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( FENDofsTest );

  CPPUNIT_TEST( testLagrangeNDofsByType );
  CPPUNIT_TEST( testMonomialNDofsByType );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testPolygonNDofs );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  typedef FE<1,LAGRANGE> Lagrange1;
  typedef FE<1,MONOMIAL> Monomial1;
#if LIBMESH_DIM > 1
  typedef FE<2,LAGRANGE> Lagrange2;
  typedef FE<2,MONOMIAL> Monomial2;
#endif
#if LIBMESH_DIM > 2
  typedef FE<3,LAGRANGE> Lagrange3;
  typedef FE<3,MONOMIAL> Monomial3;
#endif

public:
  void setUp() {}
  void tearDown() {}

  void testLagrangeNDofsByType()
  {
    LOG_UNIT_TEST;

    // A Lagrange element has a dof per node, so these are node counts
    CPPUNIT_ASSERT_EQUAL(2u, Lagrange1::n_dofs(EDGE2, FIRST));
    CPPUNIT_ASSERT_EQUAL(3u, Lagrange1::n_dofs(EDGE3, SECOND));
#if LIBMESH_DIM > 1
    CPPUNIT_ASSERT_EQUAL(3u, Lagrange2::n_dofs(TRI6, FIRST));
    CPPUNIT_ASSERT_EQUAL(6u, Lagrange2::n_dofs(TRI6, SECOND));
    CPPUNIT_ASSERT_EQUAL(8u, Lagrange2::n_dofs(QUAD8, SECOND));
#endif
#if LIBMESH_DIM > 2
    CPPUNIT_ASSERT_EQUAL(20u, Lagrange3::n_dofs(HEX20, SECOND));
    CPPUNIT_ASSERT_EQUAL(27u, Lagrange3::n_dofs(HEX27, SECOND));
#endif

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Lagrange has no fourth-order member, and a hex is not a shape a
    // constant approximation can live on
    CPPUNIT_ASSERT_THROW_MESSAGE("Lagrange FOURTH order not rejected",
                                 Lagrange1::n_dofs(EDGE2, FOURTH),
                                 libMesh::LogicError);
#if LIBMESH_DIM > 2
    CPPUNIT_ASSERT_THROW_MESSAGE("CONSTANT Lagrange on a Hex not rejected",
                                 Lagrange3::n_dofs(HEX8, CONSTANT),
                                 libMesh::LogicError);
#endif
#endif
  }

  void testMonomialNDofsByType()
  {
    LOG_UNIT_TEST;

    // A Monomial basis has as many dofs as there are monomials of the
    // given degree in the element's dimension, whatever the shape
    CPPUNIT_ASSERT_EQUAL(1u, Monomial1::n_dofs(EDGE2, CONSTANT));
    CPPUNIT_ASSERT_EQUAL(2u, Monomial1::n_dofs(EDGE2, FIRST));
#if LIBMESH_DIM > 1
    CPPUNIT_ASSERT_EQUAL(3u, Monomial2::n_dofs(TRI3, FIRST));
    CPPUNIT_ASSERT_EQUAL(6u, Monomial2::n_dofs(QUAD4, SECOND));
#endif
#if LIBMESH_DIM > 2
    // Unlike Lagrange, Monomial answers for any order: in 3D the count
    // is the number of monomials of degree o, (o+1)(o+2)(o+3)/6
    CPPUNIT_ASSERT_EQUAL(35u, Monomial3::n_dofs(HEX8, FOURTH));
    CPPUNIT_ASSERT_EQUAL(56u, Monomial3::n_dofs(HEX8, FIFTH));
#endif
  }

#if LIBMESH_DIM > 1
  void testPolygonNDofs()
  {
    LOG_UNIT_TEST;

    // A polygon's dof count is its node count, which its element type
    // does not fix, so only the overload taking an Elem can answer, and
    // only at first order.
    C0Polygon pentagon(5);

    CPPUNIT_ASSERT_EQUAL(pentagon.n_nodes(),
                         Lagrange2::n_dofs(&pentagon, FIRST));

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    CPPUNIT_ASSERT_THROW_MESSAGE("A polygon has no second-order Lagrange basis",
                                 Lagrange2::n_dofs(&pentagon, SECOND),
                                 libMesh::LogicError);
    CPPUNIT_ASSERT_THROW_MESSAGE("A polygon has no constant Lagrange basis",
                                 Lagrange2::n_dofs(&pentagon, CONSTANT),
                                 libMesh::LogicError);
    CPPUNIT_ASSERT_THROW_MESSAGE("An ElemType alone cannot count a polygon's dofs",
                                 Lagrange2::n_dofs(C0POLYGON, FIRST),
                                 libMesh::LogicError);
#endif
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( FENDofsTest );
