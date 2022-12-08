#include <libmesh/elem.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/quadrature.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/utility.h>

#include <iomanip>
#include <numeric> // std::iota

#include "libmesh_cppunit.h"


using namespace libMesh;

#define MACROCOMMA ,

#if LIBMESH_DIM > 2
#define TEST_ONE_ORDER(qtype, order, maxorder)                          \
  CPPUNIT_TEST( testBuild<qtype MACROCOMMA order> );                    \
  CPPUNIT_TEST( test1DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> ); \
  CPPUNIT_TEST( test2DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> ); \
  CPPUNIT_TEST( test3DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> );
#elif LIBMESH_DIM > 1
#define TEST_ONE_ORDER(qtype, order, maxorder)                          \
  CPPUNIT_TEST( testBuild<qtype MACROCOMMA order> );                    \
  CPPUNIT_TEST( test1DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> ); \
  CPPUNIT_TEST( test2DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> );
#else
#define TEST_ONE_ORDER(qtype, order, maxorder)                          \
  CPPUNIT_TEST( testBuild<qtype MACROCOMMA order> );                    \
  CPPUNIT_TEST( test1DWeights<qtype MACROCOMMA order MACROCOMMA maxorder> );
#endif

// std::min isn't constexpr, and C++03 lacks constexpr anyway
#define mymin(a, b) (a < b ? a : b)

#define TEST_NINTH_ORDER(qtype, maxorder)               \
  TEST_ONE_ORDER(qtype, FIRST, mymin(1,maxorder));      \
  TEST_ONE_ORDER(qtype, SECOND, mymin(2,maxorder));     \
  TEST_ONE_ORDER(qtype, THIRD, mymin(3,maxorder));      \
  TEST_ONE_ORDER(qtype, FOURTH, mymin(4,maxorder));     \
  TEST_ONE_ORDER(qtype, FIFTH, mymin(5,maxorder));      \
  TEST_ONE_ORDER(qtype, SIXTH, mymin(6,maxorder));      \
  TEST_ONE_ORDER(qtype, SEVENTH, mymin(7,maxorder));    \
  TEST_ONE_ORDER(qtype, EIGHTH, mymin(8,maxorder));     \
  TEST_ONE_ORDER(qtype, NINTH, mymin(9,maxorder));

#define TEST_TWENTIETH_ORDER(qtype, maxorder)           \
  TEST_NINTH_ORDER(qtype, maxorder)                     \
  TEST_ONE_ORDER(qtype, TENTH, mymin(10,maxorder));     \
  TEST_ONE_ORDER(qtype, ELEVENTH, mymin(11,maxorder));  \
  TEST_ONE_ORDER(qtype, TWELFTH, mymin(12,maxorder));   \
  TEST_ONE_ORDER(qtype, THIRTEENTH, mymin(13,maxorder));\
  TEST_ONE_ORDER(qtype, FOURTEENTH, mymin(14,maxorder));\
  TEST_ONE_ORDER(qtype, FIFTEENTH, mymin(15,maxorder)); \
  TEST_ONE_ORDER(qtype, SIXTEENTH, mymin(16,maxorder)); \
  TEST_ONE_ORDER(qtype, SEVENTEENTH, mymin(17,maxorder));\
  TEST_ONE_ORDER(qtype, EIGHTEENTH, mymin(18,maxorder));\
  TEST_ONE_ORDER(qtype, NINETEENTH, mymin(19,maxorder));\
  TEST_ONE_ORDER(qtype, TWENTIETH, mymin(20,maxorder));

#define LIBMESH_ASSERT_REALS_EQUAL(first, second, tolerance)            \
  if (std::abs(first-second) >= tolerance)                              \
    {                                                                   \
      std::cerr << "first = " << first << std::endl;                    \
      std::cerr << "second = " << second << std::endl;                  \
      std::cerr << "error = " << std::abs(first-second) << std::endl;   \
      std::cerr << "tolerance = " << tolerance << std::endl;            \
    }                                                                   \
  CPPUNIT_ASSERT (std::abs(first-second) < tolerance)

class QuadratureTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( QuadratureTest );

  TEST_TWENTIETH_ORDER(QGAUSS, 9999);
  TEST_ONE_ORDER(QSIMPSON, FIRST,  1);
  TEST_ONE_ORDER(QSIMPSON, SECOND, 2);
  TEST_ONE_ORDER(QSIMPSON, THIRD,  3);
  TEST_ONE_ORDER(QTRAP, FIRST, 1);
  TEST_NINTH_ORDER(QGRID, 1);

  // In general, QNodal rules (e.g. QTRAP) are only exact for linears.
  // QSIMPSON is a special case of a nodal quadrature which obtains
  // higher accuracy.
  TEST_ONE_ORDER(QNODAL, /*ignored*/FIRST, /*max order=*/1);

  TEST_NINTH_ORDER(QGAUSS_LOBATTO, 9);

  // The super-high-order Gauss-Lobatto quadrature rules only take
  // ~0.2s to test for me these days. - RHS
  TEST_ONE_ORDER(QGAUSS_LOBATTO, ELEVENTH,    11);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTEENTH,  13);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, FIFTEENTH,   15);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, SEVENTEENTH, 17);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, NINETEENTH,  19);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYFIRST, 21);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYTHIRD, 23);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYFIFTH, 25);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYSEVENTH, 27);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, TWENTYNINTH, 29);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYFIRST, 31);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYTHIRD, 33);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYFIFTH, 35);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYSEVENTH, 37);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, THIRTYNINTH, 39);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, FORTYFIRST, 41);
  TEST_ONE_ORDER(QGAUSS_LOBATTO, FORTYTHIRD, 43);

  // Might as well test super-high 1-D Gauss while we're at it.
  // Some of the really high order stuff is getting pulled in
  // indirectly; we'll just hit what lcov says we were missing.
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA TWENTYFIRST MACROCOMMA 21> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA TWENTYFIFTH MACROCOMMA 25> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA TWENTYSEVENTH MACROCOMMA 27> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA TWENTYNINTH MACROCOMMA 29> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA THIRTYFIRST MACROCOMMA 31> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA THIRTYTHIRD MACROCOMMA 33> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA THIRTYFIFTH MACROCOMMA 35> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA THIRTYSEVENTH MACROCOMMA 37> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA THIRTYNINTH MACROCOMMA 39> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA FORTYFIRST MACROCOMMA 41> );
  CPPUNIT_TEST( test1DWeights<QGAUSS MACROCOMMA FORTYTHIRD MACROCOMMA 43> );

  // Test monomial quadrature rules on quads and hexes
  CPPUNIT_TEST( testMonomialQuadrature );

  // Test nodal quadrature rules on every specific type
  CPPUNIT_TEST( testNodalQuadrature );

  // Test quadrature rules on Triangles
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testTriQuadrature );
#endif

  // Test quadrature rules on Tetrahedra
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testTetQuadrature );
#endif

  // Test Jacobi quadrature rules with special weighting function
  CPPUNIT_TEST( testJacobi );

  CPPUNIT_TEST_SUITE_END();

private:

  Real quadrature_tolerance;

  void testPolynomial(QBase & qrule,
                      int xp, int yp, int zp,
                      Real true_value)
  {
    Real sum = 0.;
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        const Point p = qrule.qp(qp);
        Real val = qrule.w(qp) * std::pow(p(0),xp);
        if (qrule.get_dim() > 1)
          val *= std::pow(p(1),yp);
        if (qrule.get_dim() > 2)
          val *= std::pow(p(2),zp);
        sum += val;
      }

    // Make sure that the weighted evaluations add up to the true
    // integral value we expect
    LIBMESH_ASSERT_REALS_EQUAL(true_value, sum, quadrature_tolerance);
  }

  void testPolynomials(QuadratureType qtype, int order,
                       ElemType elem_type,
                       const std::function<Real(int,int,int)> & true_value,
                       int exactorder)
  {
    unsigned int dim = Elem::build(elem_type)->dim();
    std::unique_ptr<QBase> qrule =
      QBase::build(qtype, dim, static_cast<Order>(order));

    // We are testing integration on the reference elements here, so
    // the element map is not relevant, and we can safely use nodal
    // Pyramid quadrature.
    if (elem_type == PYRAMID5 || elem_type == PYRAMID13 ||
        elem_type == PYRAMID14 || elem_type == PYRAMID18)
      qrule->allow_nodal_pyramid_quadrature = true;

    qrule->init (elem_type);

    const int max_y_order = dim>1 ? 0 : exactorder;
    const int max_z_order = dim>2 ? 0 : exactorder;

    for (int xp=0; xp <= exactorder; ++xp)
      for (int yp=0; yp <= max_y_order; ++yp)
        for (int zp=0; zp <= max_z_order; ++zp)
          {
            // Only try to integrate polynomials we can integrate exactly
            if (xp + yp + zp > exactorder)
              continue;

            const Real exact = true_value(xp, yp, zp);
            testPolynomial(*qrule, xp, yp, zp, exact);
          }
  }

  const std::function<Real(int,int,int)> edge_integrals =
  [](int mode, int, int) {
    return (mode % 2) ?  0 : (Real(2.0) / (mode+1));
  };

  const std::function<Real(int,int,int)> quad_integrals =
  [](int modex, int modey, int) {
    const Real exactx = (modex % 2) ?
      0 : (Real(2.0) / (modex+1));

    const Real exacty = (modey % 2) ?
      0 : (Real(2.0) / (modey+1));

    return exactx*exacty;
  };

  const std::function<Real(int,int,int)> tri_integrals =
  [](int x_power, int y_power, int) {
    // Compute the true integral, a! b! / (a + b + 2)!
    Real analytical = 1.0;

    unsigned
      larger_power = std::max(x_power, y_power),
      smaller_power = std::min(x_power, y_power);

    // Cancel the larger of the two numerator terms with the
    // denominator, and fill in the remaining entries.
    std::vector<unsigned>
      numerator(smaller_power > 1 ? smaller_power-1 : 0),
      denominator(2+smaller_power);

    // Fill up the vectors with sequences starting at the right values.
    std::iota(numerator.begin(), numerator.end(), 2);
    std::iota(denominator.begin(), denominator.end(), larger_power+1);

    // The denominator is guaranteed to have more terms...
    for (std::size_t i=0; i<denominator.size(); ++i)
      {
        if (i < numerator.size())
          analytical *= numerator[i];
        analytical /= denominator[i];
      }
    return analytical;
  };

  const std::function<Real(int,int,int)> hex_integrals =
  [](int modex, int modey, int modez) {
    const Real exactx = (modex % 2) ?
      0 : (Real(2.0) / (modex+1));

    const Real exacty = (modey % 2) ?
      0 : (Real(2.0) / (modey+1));

    const Real exactz = (modez % 2) ?
      0 : (Real(2.0) / (modez+1));

    return exactx*exacty*exactz;
  };

  const std::function<Real(int,int,int)> tet_integrals =
  [](int x_power, int y_power, int z_power) {
    // Compute the true integral, a! b! c! / (a + b + c + 3)!
    Real analytical = 1.0;

    // Sort the a, b, c values
    int sorted_powers[3] = {x_power, y_power, z_power};
    std::sort(sorted_powers, sorted_powers+3);

    // Cancel the largest power with the denominator, fill in the
    // entries for the remaining numerator terms and the denominator.
    std::vector<int>
      numerator_1(sorted_powers[0] > 1 ? sorted_powers[0]-1 : 0),
      numerator_2(sorted_powers[1] > 1 ? sorted_powers[1]-1 : 0),
      denominator(3 + sorted_powers[0] + sorted_powers[1]);

    // Fill up the vectors with sequences starting at the right values.
    std::iota(numerator_1.begin(), numerator_1.end(), 2);
    std::iota(numerator_2.begin(), numerator_2.end(), 2);
    std::iota(denominator.begin(), denominator.end(), sorted_powers[2]+1);

    // The denominator is guaranteed to have the most terms...
    for (std::size_t i=0; i<denominator.size(); ++i)
      {
        if (i < numerator_1.size())
          analytical *= numerator_1[i];

        if (i < numerator_2.size())
          analytical *= numerator_2[i];

        analytical /= denominator[i];
      }
    return analytical;
  };

  const std::function<Real(int,int,int)> prism_integrals =
  [this](int modex, int modey, int modez) {
    const Real exactz = (modez % 2) ?
      0 : (Real(2.0) / (modez+1));

    return exactz * tri_integrals(modex, modey, 0);
  };

  const std::function<Real(int,int,int)> pyramid_integrals =
  [](int modex, int modey, int modez) {

    const int binom = Utility::binomial(modex+modey+modez+3, modez);

    if (modex%2 || modey%2)
      return Real(0);

    return Real(4)/((modex+1)*(modey+1)*binom*(modex+modey+modez+3));
  };


public:
  void setUp ()
  { quadrature_tolerance = TOLERANCE * std::sqrt(TOLERANCE); }

  void tearDown ()
  {}

  void testNodalQuadrature ()
  {
    LOG_UNIT_TEST;

    const std::vector<std::vector<ElemType>> all_types =
      {{EDGE2, EDGE3, EDGE4},
       {TRI3, TRI6, TRI7},
       {QUAD4, QUAD8, QUAD9},
       {TET4, TET10, TET14},
       {HEX8, HEX20, HEX27},
       {PRISM6, PRISM15, PRISM18, PRISM20, PRISM21},
       {PYRAMID5, PYRAMID13, PYRAMID14, PYRAMID18}};

    const std::function<Real(int,int,int)> true_values[] =
      {edge_integrals,
       tri_integrals,
       quad_integrals,
       tet_integrals,
       hex_integrals,
       prism_integrals,
       pyramid_integrals};

    for (auto i : index_range(all_types))
      for (ElemType elem_type : all_types[i])
        {
          const unsigned int order = Elem::build(elem_type)->default_order();

          unsigned int exactorder = order;
          // There are some sad exceptions to the "positive nodal
          // quadrature can integrate polynomials up to the element
          // order" rule.
          //
          // Some "quadratic" elements can only exactly integrate
          // linears when mass lumping.
          if (elem_type == TRI6 || elem_type == QUAD8 ||
              elem_type == TET10 || elem_type == HEX20 ||
              elem_type == PRISM15 || elem_type == PRISM18 ||
              elem_type == PYRAMID13 || elem_type == PYRAMID14 ||
          // And some partially-cubic elements can only do quadratics
              elem_type == TET14 || elem_type == PRISM20)
            exactorder--;

          // And one partially-cubic element can only do linears.
          // Seriously.
          if (elem_type == PYRAMID18)
            exactorder=1;

          testPolynomials(QNODAL, order, elem_type, true_values[i], exactorder);
        }
  }


  void testMonomialQuadrature ()
  {
    LOG_UNIT_TEST;

    const std::vector<std::vector<ElemType>> all_types =
      {{EDGE2}, {QUAD4, TRI3}, {HEX8, TET4, PRISM6, PYRAMID5}};
    const std::vector<std::vector<std::function<Real(int,int,int)>>>
      true_values =
      {{edge_integrals},
       {quad_integrals, tri_integrals},
       {hex_integrals, tet_integrals, prism_integrals, pyramid_integrals}};

    for (auto i : index_range(all_types))
      for (auto j : index_range(all_types[i]))
        for (int order=0; order<17; ++order)
          testPolynomials(QMONOMIAL, order, all_types[i][j], true_values[i][j], order);
  }

  void testTetQuadrature ()
  {
    LOG_UNIT_TEST;

    // There are 3 different families of quadrature rules for tetrahedra
    QuadratureType qtype[3] = {QCONICAL, QGRUNDMANN_MOLLER, QGAUSS};

    int end_order = 7;
    // Our higher order tet rules were only computed to double precision
    if (quadrature_tolerance < 1e-16)
      end_order = 2;

    for (int qt=0; qt<3; ++qt)
      for (int order=0; order<end_order; ++order)
        testPolynomials(qtype[qt], order, TET4, tet_integrals, order);
  }

  void testTriQuadrature ()
  {
    LOG_UNIT_TEST;

    QuadratureType qtype[4] = {QCONICAL, QCLOUGH, QGAUSS, QGRUNDMANN_MOLLER};

    for (int qt=0; qt<4; ++qt)
      for (int order=0; order<10; ++order)
        testPolynomials(qtype[qt], order, TRI3, tri_integrals, order);
  }

  void testJacobi ()
  {
    LOG_UNIT_TEST;

    constexpr int max_order = 43;  // We seriously implemented that?!

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
    std::vector<std::vector<Real>> true_integrals(2);

    // alpha=1 integral values
    // int((1-x)*x^p, x=0..1) = 1 / (p^2 + 3p + 2)
    true_integrals[0].resize(max_order);
    for (std::size_t p=0; p<true_integrals[0].size(); ++p)
      true_integrals[0][p] = 1. / (p*p + 3.*p + 2.);

    // alpha=2 integral values
    // int((1-x)^2*x^p, x=0..1) = 2 / (p^3 + 6*p^2 + 11*p + 6)
    true_integrals[1].resize(max_order);
    for (std::size_t p=0; p<true_integrals[1].size(); ++p)
      true_integrals[1][p] = 2. / (p*p*p + 6.*p*p + 11.*p + 6.);

    // Test both types of Jacobi quadrature rules
    for (int qt=0; qt<2; ++qt)
      {
        for (int order=0; order<max_order; ++order)
          {
            std::unique_ptr<QBase> qrule = QBase::build(qtype[qt],
                                                        /*dim=*/1,
                                                        static_cast<Order>(order));

            // Initialize on a 1D element, EDGE2/3/4 should not matter...
            qrule->init (EDGE2);

            // Test the sum of the weights for this order
            Real sumw = 0.;
            for (unsigned int qp=0; qp<qrule->n_points(); qp++)
              sumw += qrule->w(qp);

            // Make sure that the weights add up to the value we expect
            LIBMESH_ASSERT_REALS_EQUAL(initial_sum_weights[qt], sumw, quadrature_tolerance);

            // Test integrating different polynomial powers
            for (int testpower=0; testpower<=order; ++testpower)
              {
                // Note that the weighting function, (1-x)^alpha *
                // (1+x)^beta, is built into these quadrature rules;
                // the polynomials we actually integrate are just the
                // usual monomials.
                Real sumq = 0.;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                  sumq += qrule->w(qp) * std::pow(qrule->qp(qp)(0), testpower);

                // Make sure that the computed integral agrees with the "true" value
                LIBMESH_ASSERT_REALS_EQUAL(true_integrals[qt][testpower], sumq, quadrature_tolerance);
              } // end for(testpower)
          } // end for(order)
      } // end for(qt)
  } // testJacobi



  template <QuadratureType qtype, Order order>
  void testBuild ()
  {
    LOG_UNIT_TEST;

    std::unique_ptr<QBase> qrule1D = QBase::build (qtype, 1, order);
    std::unique_ptr<QBase> qrule2D = QBase::build (qtype, 2, order);
    std::unique_ptr<QBase> qrule3D = QBase::build (qtype, 3, order);

    CPPUNIT_ASSERT_EQUAL ( qtype, qrule1D->type() );
    CPPUNIT_ASSERT_EQUAL ( qtype, qrule2D->type() );
    CPPUNIT_ASSERT_EQUAL ( qtype, qrule3D->type() );

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
    LOG_UNIT_TEST;

    testPolynomials(qtype, order, EDGE3, edge_integrals, exactorder);
  }



  //-------------------------------------------------------
  // 2D Quadrature Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test2DWeights ()
  {
    LOG_UNIT_TEST;

    testPolynomials(qtype, order, QUAD8, quad_integrals, exactorder);

    // We may eventually support Gauss-Lobatto type quadrature on triangles...
    if (qtype == QGAUSS_LOBATTO)
      return;

    // QGrid needs to be changed to use symmetric offsets on triangles
    // so it can at *least* get linears right...
    if (qtype == QGRID)
      testPolynomials(qtype, order, TRI6, tri_integrals, 0);
    // QSimpson doesn't even get all quadratics on a triangle
    else if (qtype == QSIMPSON)
      testPolynomials(qtype, order, TRI6, tri_integrals, std::min(1u,exactorder));
    else
      testPolynomials(qtype, order, TRI6, tri_integrals, exactorder);
  }



  //-------------------------------------------------------
  // 3D Gauss Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test3DWeights ()
  {
    LOG_UNIT_TEST;

    testPolynomials(qtype, order, HEX20, hex_integrals, exactorder);

    // We may eventually support Gauss-Lobatto type quadrature on tets and prisms...
    if (qtype == QGAUSS_LOBATTO)
      return;

    // QGrid needs to be changed to use symmetric offsets on
    // non-tensor product elements so it can at *least* get linears
    // right...
    if (qtype == QGRID)
      {
        testPolynomials(qtype, order, TET10, tet_integrals, 0);
        testPolynomials(qtype, order, PYRAMID14, pyramid_integrals, 0);
      }
    // QSimpson doesn't get all quadratics on a simplex or its extrusion
    else if (qtype == QSIMPSON)
      {
        testPolynomials(qtype, order, TET10, tet_integrals, std::min(1u,exactorder));
        testPolynomials(qtype, order, PRISM15, prism_integrals, std::min(1u,exactorder));

        // And on pyramids we gave up and redid QTrap
        testPolynomials(qtype, order, PYRAMID14, pyramid_integrals, std::min(1u,exactorder));
      }
    else
      {
        testPolynomials(qtype, order, TET10, tet_integrals, exactorder);
        testPolynomials(qtype, order, PRISM15, prism_integrals, exactorder);
        testPolynomials(qtype, order, PYRAMID14, pyramid_integrals, exactorder);
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( QuadratureTest );
