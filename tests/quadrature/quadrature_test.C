#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <libmesh/quadrature.h>

using namespace libMesh;

#define MACROCOMMA ,

#define TEST_ONE_ORDER(qtype, order) \
  CPPUNIT_TEST( testBuild<qtype MACROCOMMA order> ); \
  CPPUNIT_TEST( test1DWeights<qtype MACROCOMMA order> ); \
  CPPUNIT_TEST( test2DWeights<qtype MACROCOMMA order> ); \
  CPPUNIT_TEST( test3DWeights<qtype MACROCOMMA order> );

#define TEST_ALL_ORDERS(qtype) \
  TEST_ONE_ORDER(qtype, FIRST); \
  TEST_ONE_ORDER(qtype, SECOND); \
  TEST_ONE_ORDER(qtype, THIRD); \
  TEST_ONE_ORDER(qtype, FOURTH); \
  TEST_ONE_ORDER(qtype, FIFTH); \
  TEST_ONE_ORDER(qtype, SIXTH); \
  TEST_ONE_ORDER(qtype, SEVENTH); \


class QuadratureTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( QuadratureTest );

  TEST_ALL_ORDERS(QGAUSS);
  TEST_ALL_ORDERS(QSIMPSON);
  TEST_ALL_ORDERS(QTRAP);
  TEST_ALL_ORDERS(QGRID);
  TEST_ALL_ORDERS(QGRUNDMANN_MOLLER);
  TEST_ALL_ORDERS(QMONOMIAL);
  TEST_ALL_ORDERS(QCONICAL);
  TEST_ALL_ORDERS(QCLOUGH);

  CPPUNIT_TEST_SUITE_END();

private:


public:
  void setUp ()
  {}

  void tearDown ()
  {}



  template <QuadratureType qtype, Order order>
  void testBuild ()
  {
    AutoPtr<QBase> qrule1D = QBase::build (qtype, 1, order);
    AutoPtr<QBase> qrule2D = QBase::build (qtype, 2, order);
    AutoPtr<QBase> qrule3D = QBase::build (qtype, 3, order);

    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(1) , qrule1D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(2) , qrule2D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(3) , qrule3D->get_dim() );
  }



  //-------------------------------------------------------
  // 1D Quadrature Rule Test
  template <QuadratureType qtype, Order order>
  void test1DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype , 1, order);
    qrule->init (EDGE3);

    Real sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , sum , TOLERANCE*TOLERANCE );
  }



  //-------------------------------------------------------
  // 2D Quadrature Rule Test
  template <QuadratureType qtype, Order order>
  void test2DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype, 2, order);
    qrule->init (QUAD8);

    Real sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0 , sum , TOLERANCE*TOLERANCE );

    qrule->init (TRI6);

    sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5 , sum , TOLERANCE*TOLERANCE );
  }



  //-------------------------------------------------------
  // 3D Gauss Rule Test
  template <QuadratureType qtype, Order order>
  void test3DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype, 3, order);
    qrule->init (HEX20);

    Real sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0 , sum , TOLERANCE*TOLERANCE );

    qrule->init (TET10);

    sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1./6., sum , TOLERANCE*TOLERANCE );

    qrule->init (PRISM15);

    sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., sum , TOLERANCE*TOLERANCE );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( QuadratureTest );
