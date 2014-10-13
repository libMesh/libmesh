// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/quadrature.h>
#include <libmesh/string_to_enum.h>

#include <iomanip>

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

  // The TEST_ALL_ORDERS macro only goes up to 9th-order
  TEST_ALL_ORDERS(QGAUSS_LOBATTO, 9);

  // The Gauss-Lobatto quadrature rules passed all these tests during
  // development, but we don't run them with every 'make check'
  // because it takes too long.
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, ELEVENTH,    11);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTEENTH,  13);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, FIFTEENTH,   15);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, SEVENTEENTH, 17);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, NINETEENTH,  19);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYFIRST, 21);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYTHIRD, 23);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYFIFTH, 25);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYSEVENTH, 27);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYNINTH, 29);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYFIRST, 31);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYTHIRD, 33);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYFIFTH, 35);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYSEVENTH, 37);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYNINTH, 39);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, FORTYFIRST, 41);
  // TEST_ONE_ORDER(QGAUSS_LOBATTO, FORTYTHIRD, 43);

// Edges/Tris only
//  TEST_ALL_ORDERS(QCLOUGH, 9999);

// Quads/Hexes only
//  TEST_ALL_ORDERS(QMONOMIAL, 1); // need "non-tensor" option?

// Tets only
//  TEST_ALL_ORDERS(QGRUNDMANN_MOLLER, 9999);

// Tris+Tets only
//  TEST_ALL_ORDERS(QCONICAL, 9999);

// Test Jacobi quadrature rules with special weighting function
  CPPUNIT_TEST( testJacobi );

  CPPUNIT_TEST_SUITE_END();

private:


public:
  void setUp ()
  {}

  void tearDown ()
  {}

  void testJacobi ()
  {
    // LibMesh supports two different types of Jacobi quadrature
    QuadratureType qtype[2] = {QJACOBI_1_0, QJACOBI_2_0};

    // The weights of the Jacobi quadrature rules in libmesh have been
    // scaled based on their intended use:
    // (alpha=1, beta=0) rule weights sum to 1/2.
    // (alpha=2, beta=0) rule weights sum to 1/3.
    Real initial_sum_weights[2] = {.5, 1./3.};

    // The points of the Jacobi rules in LibMesh are also defined on
    // [0,1]... this has to be taken into account when computing the
    // exact integral values in Maple!  Also note: if you scale the
    // points to a different interval, you need to also compute what
    // the sum of the weights should be for that interval, it will not
    // simply be the element length for weighted quadrature rules like
    // Jacobi.  For general alpha and beta=0, the sum of the weights
    // on the interval [-1,1] is 2^(alpha+1) / (alpha+1).
    std::vector<std::vector<Real> > true_integrals(2);

    // alpha=1 integral values
    // int((1-x)*x^p, x=0..1) = 1 / (p^2 + 3p + 2)
    true_integrals[0].resize(10);
    for (unsigned p=0; p<true_integrals[0].size(); ++p)
      true_integrals[0][p] = 1. / (p*p + 3.*p + 2.);

    // alpha=2 integral values
    // int((1-x)^2*x^p, x=0..1) = 2 / (p^3 + 6*p^2 + 11*p + 6)
    true_integrals[1].resize(10);
    for (unsigned p=0; p<true_integrals[1].size(); ++p)
      true_integrals[1][p] = 2. / (p*p*p + 6.*p*p + 11.*p + 6.);

    // Test both types of Jacobi quadrature rules
    for (int qt=0; qt<2; ++qt)
      {
        for (int order=0; order<10; ++order)
          {
            AutoPtr<QBase> qrule = QBase::build(qtype[qt],
                                                /*dim=*/1,
                                                static_cast<Order>(order));

            // Initialize on a 1D element, EDGE2/3/4 should not matter...
            qrule->init (EDGE2);

            // Test the sum of the weights for this order
            Real sumw = 0.;
            for (unsigned int qp=0; qp<qrule->n_points(); qp++)
              sumw += qrule->w(qp);

            // Make sure that the weights add up to the value we expect
            CPPUNIT_ASSERT_DOUBLES_EQUAL(initial_sum_weights[qt], sumw, TOLERANCE*TOLERANCE);

            // Test integrating different polynomial powers
            for (int testpower=0; testpower<=order; ++testpower)
              {
                // Don't try testpowers larger than the order of the method
                if (testpower > order)
                  continue;

                // Note that the weighting function, (1-x)^alpha *
                // (1+x)^beta, is built into these quadrature rules;
                // the polynomials we actually integrate are just the
                // usual monomials.
                Real sumq = 0.;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                  sumq += qrule->w(qp) * std::pow(qrule->qp(qp)(0), testpower);

                // Make sure that the computed integral agrees with the "true" value
                CPPUNIT_ASSERT_DOUBLES_EQUAL(true_integrals[qt][testpower], sumq, TOLERANCE*TOLERANCE);
              } // end for(testpower)
          } // end for(order)
      } // end for(qt)
  } // testJacobi



  template <QuadratureType qtype, Order order>
  void testBuild ()
  {
    AutoPtr<QBase> qrule1D = QBase::build (qtype, 1, order);
    AutoPtr<QBase> qrule2D = QBase::build (qtype, 2, order);
    AutoPtr<QBase> qrule3D = QBase::build (qtype, 3, order);

    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(1) , qrule1D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(2) , qrule2D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(3) , qrule3D->get_dim() );

    // Test the enum-to-string conversion for this qtype is
    // implemented, but not what the actual value is.
    Utility::enum_to_string(qtype);
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

    // We may eventually support Gauss-Lobatto type quadrature on triangles...
    if (qtype != QGAUSS_LOBATTO)
      {
        qrule->init (TRI6);

        Real sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp);

        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5 , sum , TOLERANCE*TOLERANCE );
      }
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

    // We may eventually support Gauss-Lobatto type quadrature on tets and prisms...
    if (qtype != QGAUSS_LOBATTO)
      {
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
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( QuadratureTest );
