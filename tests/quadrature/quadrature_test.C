// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/quadrature.h>

using namespace libMesh;

#define MACROCOMMA ,

#define TEST_ONE_ORDER(qtype, order, maxorder) \
  CPPUNIT_TEST( testBuild<qtype MACROCOMMA order> ); \
  CPPUNIT_TEST( test1DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> ); \
  CPPUNIT_TEST( test2DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> ); \
  CPPUNIT_TEST( test3DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> );

// std::min isn't constexpr, and C++03 lacks constexpr anyway
#define mymin(a, b) (a < b ? a : b)

#define TEST_ALL_ORDERS(qtype, maxorder) \
  TEST_ONE_ORDER(qtype, FIRST, mymin(1,maxorder)); \
  TEST_ONE_ORDER(qtype, SECOND, mymin(2,maxorder)); \
  TEST_ONE_ORDER(qtype, THIRD, mymin(3,maxorder)); \
  TEST_ONE_ORDER(qtype, FOURTH, mymin(4,maxorder)); \
  TEST_ONE_ORDER(qtype, FIFTH, mymin(5,maxorder)); \
  TEST_ONE_ORDER(qtype, SIXTH, mymin(6,maxorder)); \
  TEST_ONE_ORDER(qtype, SEVENTH, mymin(7,maxorder)); \
  TEST_ONE_ORDER(qtype, EIGHTH, mymin(8,maxorder)); \
  TEST_ONE_ORDER(qtype, NINTH, mymin(9,maxorder));


class QuadratureTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( QuadratureTest );

  TEST_ALL_ORDERS(QGAUSS, 9999);
  TEST_ONE_ORDER(QSIMPSON, FIRST,  3);
  TEST_ONE_ORDER(QSIMPSON, SECOND, 3);
  TEST_ONE_ORDER(QSIMPSON, THIRD,  3);
  TEST_ONE_ORDER(QTRAP, FIRST, 1);
  TEST_ALL_ORDERS(QGRID, 1);

// Edges/Tris only
//  TEST_ALL_ORDERS(QCLOUGH, 9999);

// Quads/Hexes only
//  TEST_ALL_ORDERS(QMONOMIAL, 1); // need "non-tensor" option?

// Tets only
//  TEST_ALL_ORDERS(QGRUNDMANN_MOLLER, 9999);

// Tris+Tets only
//  TEST_ALL_ORDERS(QCONICAL, 9999);

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
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test1DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype , 1, order);
    qrule->init (EDGE3);

    for (unsigned int mode=0; mode <= exactorder; ++mode)
      {
        Real sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp) * std::pow(qrule->qp(qp)(0), static_cast<Real>(mode));

        const Real exact = (mode % 2) ?
          0 : (Real(2.0) / (mode+1));

if (std::abs(exact - sum) >= TOLERANCE*TOLERANCE)
{
    std::cout << "qtype = " << qtype << std::endl;
    std::cout << "order = " << order << std::endl;
    std::cout << "exactorder = " << exactorder << std::endl;
    std::cout << "mode = " << mode << std::endl;
    std::cout << "exact = " << exact << std::endl;
    std::cout << "sum = " << sum << std::endl << std::endl;
}

        CPPUNIT_ASSERT_DOUBLES_EQUAL( exact , sum , TOLERANCE*TOLERANCE );
      }
  }



  //-------------------------------------------------------
  // 2D Quadrature Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test2DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype, 2, order);
    qrule->init (QUAD8);

    for (unsigned int modex=0; modex <= exactorder; ++modex)
      for (unsigned int modey=0; modey <= exactorder; ++modey)
        {
          Real sum = 0;

          for (unsigned int qp=0; qp<qrule->n_points(); qp++)
            sum += qrule->w(qp) * std::pow(qrule->qp(qp)(0), static_cast<Real>(modex))
                                * std::pow(qrule->qp(qp)(1), static_cast<Real>(modey));

          const Real exactx = (modex % 2) ?
            0 : (Real(2.0) / (modex+1));

          const Real exacty = (modey % 2) ?
            0 : (Real(2.0) / (modey+1));

          const Real exact = exactx*exacty;

          CPPUNIT_ASSERT_DOUBLES_EQUAL( exact , sum , TOLERANCE*TOLERANCE );
        }

    qrule->init (TRI6);

    Real sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5 , sum , TOLERANCE*TOLERANCE );
  }



  //-------------------------------------------------------
  // 3D Gauss Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test3DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(qtype, 3, order);
    qrule->init (HEX20);

    for (unsigned int modex=0; modex <= exactorder; ++modex)
      for (unsigned int modey=0; modey <= exactorder; ++modey)
        for (unsigned int modez=0; modez <= exactorder; ++modez)
          {
            Real sum = 0;

            for (unsigned int qp=0; qp<qrule->n_points(); qp++)
              sum += qrule->w(qp) * std::pow(qrule->qp(qp)(0), static_cast<Real>(modex))
                                  * std::pow(qrule->qp(qp)(1), static_cast<Real>(modey))
                                  * std::pow(qrule->qp(qp)(2), static_cast<Real>(modez));

            const Real exactx = (modex % 2) ?
              0 : (Real(2.0) / (modex+1));

            const Real exacty = (modey % 2) ?
              0 : (Real(2.0) / (modey+1));

            const Real exactz = (modez % 2) ?
              0 : (Real(2.0) / (modez+1));

            const Real exact = exactx*exacty*exactz;

            CPPUNIT_ASSERT_DOUBLES_EQUAL( exact , sum , TOLERANCE*TOLERANCE );
          }

    qrule->init (TET10);

    Real sum = 0;

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
