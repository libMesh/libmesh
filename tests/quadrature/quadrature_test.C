// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/quadrature.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/utility.h>

#include <iomanip>

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

#define LIBMESH_ASSERT_REALS_EQUAL(first, second, tolerance) \
  if (std::abs(first-second) >= tolerance) \
    { \
      std::cerr << "first = " << first << std::endl; \
      std::cerr << "second = " << second << std::endl; \
      std::cerr << "error = " << std::abs(first-second) << std::endl; \
      std::cerr << "tolerance = " << tolerance << std::endl; \
    } \
  CPPUNIT_ASSERT (std::abs(first-second) < tolerance)

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

  // Test monomial quadrature rules on quads and hexes
  CPPUNIT_TEST( testMonomialQuadrature );

  // Test quadrature rules on Triangles
  CPPUNIT_TEST( testTriQuadrature );

  // Test quadrature rules on Tetrahedra
  CPPUNIT_TEST( testTetQuadrature );

  // Test Jacobi quadrature rules with special weighting function
  CPPUNIT_TEST( testJacobi );

  CPPUNIT_TEST_SUITE_END();

private:

  Real quadrature_tolerance;


public:
  void setUp ()
  { quadrature_tolerance = TOLERANCE * std::sqrt(TOLERANCE); }

  void tearDown ()
  {}

  void testMonomialQuadrature ()
  {
    ElemType elem_type[2] = {QUAD4, HEX8};
    int dims[2]           = {2, 3};

    for (int i=0; i<2; ++i)
      {
        // std::cout << "\nChecking monomial quadrature on element type " << elem_type[i] << std::endl;

        for (int order=0; order<7; ++order)
          {
            UniquePtr<QBase> qrule = QBase::build(QMONOMIAL,
                                                  dims[i],
                                                  static_cast<Order>(order));
            qrule->init(elem_type[i]);

            // In 3D, max(z_power)==order, in 2D max(z_power)==0
            int max_z_power = dims[i]==2 ? 0 : order;

            for (int x_power=0; x_power<=order; ++x_power)
              for (int y_power=0; y_power<=order; ++y_power)
                for (int z_power=0; z_power<=max_z_power; ++z_power)
                  {
                    // Only try to integrate polynomials we can integrate exactly
                    if (x_power + y_power + z_power > order)
                      continue;

                    // Compute the integral via quadrature.  Note that
                    // std::pow(0,0) returns 1 in the 2D case.
                    Real sumq = 0.;
                    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                      sumq += qrule->w(qp)
                        * std::pow(qrule->qp(qp)(0), x_power)
                        * std::pow(qrule->qp(qp)(1), y_power)
                        * std::pow(qrule->qp(qp)(2), z_power);

                    // std::cout << "Quadrature of x^" << x_power
                    //           << " y^" << y_power
                    //           << " z^" << z_power
                    //           << " = " << sumq << std::endl;

                    // Copy-pasted code from test3DWeights()
                    Real exact_x = (x_power % 2) ? 0 : (Real(2.0) / (x_power+1));
                    Real exact_y = (y_power % 2) ? 0 : (Real(2.0) / (y_power+1));
                    Real exact_z = (z_power % 2) ? 0 : (Real(2.0) / (z_power+1));

                    // Handle 2D
                    if (dims[i]==2)
                      exact_z = 1.0;

                    Real exact = exact_x*exact_y*exact_z;

                    // Make sure that the quadrature solution matches the exact solution
                    LIBMESH_ASSERT_REALS_EQUAL(exact, sumq, quadrature_tolerance);
                  }
          } // end for (order)
      } // end for (i)
  }

  void testTetQuadrature ()
  {
    // There are 3 different families of quadrature rules for tetrahedra
    QuadratureType qtype[3] = {QCONICAL, QGRUNDMANN_MOLLER, QGAUSS};

    for (int qt=0; qt<3; ++qt)
      for (int order=0; order<7; ++order)
        {
          UniquePtr<QBase> qrule = QBase::build(qtype[qt],
                                                /*dim=*/3,
                                                static_cast<Order>(order));

          // Initialize on a TET element
          qrule->init (TET4);

          // Test the sum of the weights for this order
          Real sumw = 0.;
          for (unsigned int qp=0; qp<qrule->n_points(); qp++)
            sumw += qrule->w(qp);

          // Make sure that the weights add up to the value we expect
          LIBMESH_ASSERT_REALS_EQUAL(1./6., sumw, quadrature_tolerance);

          // Test integrating different polynomial powers
          for (int x_power=0; x_power<=order; ++x_power)
            for (int y_power=0; y_power<=order; ++y_power)
              for (int z_power=0; z_power<=order; ++z_power)
                {
                  // Only try to integrate polynomials we can integrate exactly
                  if (x_power + y_power + z_power > order)
                    continue;

                  // Compute the integral via quadrature
                  Real sumq = 0.;
                  for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                    sumq += qrule->w(qp)
                      * std::pow(qrule->qp(qp)(0), x_power)
                      * std::pow(qrule->qp(qp)(1), y_power)
                      * std::pow(qrule->qp(qp)(2), z_power);

                  // std::cout << "sumq = " << sumq << std::endl;

                  // Compute the true integral, a! b! c! / (a + b + c + 3)!
                  Real analytical = 1.0;
                  {
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
                    Utility::iota(numerator_1.begin(), numerator_1.end(), 2);
                    Utility::iota(numerator_2.begin(), numerator_2.end(), 2);
                    Utility::iota(denominator.begin(), denominator.end(), sorted_powers[2]+1);

                    // The denominator is guaranteed to have the most terms...
                    for (std::size_t i=0; i<denominator.size(); ++i)
                      {
                        if (i < numerator_1.size())
                          analytical *= numerator_1[i];

                        if (i < numerator_2.size())
                          analytical *= numerator_2[i];

                        analytical /= denominator[i];
                      }
                  }

                  // std::cout << "analytical = " << analytical << std::endl;

                  // Make sure that the computed integral agrees with the "true" value
                  LIBMESH_ASSERT_REALS_EQUAL(analytical, sumq, quadrature_tolerance);
                } // end for(testpower)
        } // end for(order)
  }

  void testTriQuadrature ()
  {
    QuadratureType qtype[4] = {QCONICAL, QCLOUGH, QGAUSS, QGRUNDMANN_MOLLER};

    for (int qt=0; qt<4; ++qt)
      for (int order=0; order<10; ++order)
        {
          UniquePtr<QBase> qrule = QBase::build(qtype[qt],
                                                /*dim=*/2,
                                                static_cast<Order>(order));

          // Initialize on a TRI element
          qrule->init (TRI3);

          // Test the sum of the weights for this order
          Real sumw = 0.;
          for (unsigned int qp=0; qp<qrule->n_points(); qp++)
            sumw += qrule->w(qp);

          // Make sure that the weights add up to the value we expect
          LIBMESH_ASSERT_REALS_EQUAL(0.5, sumw, quadrature_tolerance);

          // Test integrating different polynomial powers
          for (int x_power=0; x_power<=order; ++x_power)
            for (int y_power=0; y_power<=order; ++y_power)
              {
                // Only try to integrate polynomials we can integrate exactly
                if (x_power + y_power > order)
                  continue;

                // Compute the integral via quadrature
                Real sumq = 0.;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                  sumq += qrule->w(qp) * std::pow(qrule->qp(qp)(0), x_power) * std::pow(qrule->qp(qp)(1), y_power);

                // std::cout << "sumq = " << sumq << std::endl;

                // Compute the true integral, a! b! / (a + b + 2)!
                Real analytical = 1.0;
                {
                  unsigned
                    larger_power = std::max(x_power, y_power),
                    smaller_power = std::min(x_power, y_power);

                  // Cancel the larger of the two numerator terms with the
                  // denominator, and fill in the remaining entries.
                  std::vector<unsigned>
                    numerator(smaller_power > 1 ? smaller_power-1 : 0),
                    denominator(2+smaller_power);

                  // Fill up the vectors with sequences starting at the right values.
                  Utility::iota(numerator.begin(), numerator.end(), 2);
                  Utility::iota(denominator.begin(), denominator.end(), larger_power+1);

                  // The denominator is guaranteed to have more terms...
                  for (std::size_t i=0; i<denominator.size(); ++i)
                    {
                      if (i < numerator.size())
                        analytical *= numerator[i];
                      analytical /= denominator[i];
                    }
                }

                // std::cout << "analytical = " << analytical << std::endl;

                // Make sure that the computed integral agrees with the "true" value
                LIBMESH_ASSERT_REALS_EQUAL(analytical, sumq, quadrature_tolerance);
              } // end for(testpower)
        } // end for(order)
  }

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
    for (std::size_t p=0; p<true_integrals[0].size(); ++p)
      true_integrals[0][p] = 1. / (p*p + 3.*p + 2.);

    // alpha=2 integral values
    // int((1-x)^2*x^p, x=0..1) = 2 / (p^3 + 6*p^2 + 11*p + 6)
    true_integrals[1].resize(10);
    for (std::size_t p=0; p<true_integrals[1].size(); ++p)
      true_integrals[1][p] = 2. / (p*p*p + 6.*p*p + 11.*p + 6.);

    // Test both types of Jacobi quadrature rules
    for (int qt=0; qt<2; ++qt)
      {
        for (int order=0; order<10; ++order)
          {
            UniquePtr<QBase> qrule = QBase::build(qtype[qt],
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
    UniquePtr<QBase> qrule1D = QBase::build (qtype, 1, order);
    UniquePtr<QBase> qrule2D = QBase::build (qtype, 2, order);
    UniquePtr<QBase> qrule3D = QBase::build (qtype, 3, order);

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
    UniquePtr<QBase> qrule = QBase::build(qtype , 1, order);
    qrule->init (EDGE3);

    for (unsigned int mode=0; mode <= exactorder; ++mode)
      {
        Real sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp) * std::pow(qrule->qp(qp)(0), static_cast<Real>(mode));

        const Real exact = (mode % 2) ?
          0 : (Real(2.0) / (mode+1));

        if (std::abs(exact - sum) >= quadrature_tolerance)
          {
            std::cout << "qtype = " << qtype << std::endl;
            std::cout << "order = " << order << std::endl;
            std::cout << "exactorder = " << exactorder << std::endl;
            std::cout << "mode = " << mode << std::endl;
            std::cout << "exact = " << exact << std::endl;
            std::cout << "sum = " << sum << std::endl << std::endl;
          }

        LIBMESH_ASSERT_REALS_EQUAL( exact , sum , quadrature_tolerance );
      }
  }



  //-------------------------------------------------------
  // 2D Quadrature Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test2DWeights ()
  {
    UniquePtr<QBase> qrule = QBase::build(qtype, 2, order);
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

          LIBMESH_ASSERT_REALS_EQUAL( exact , sum , quadrature_tolerance );
        }

    // We may eventually support Gauss-Lobatto type quadrature on triangles...
    if (qtype != QGAUSS_LOBATTO)
      {
        qrule->init (TRI6);

        Real sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp);

        LIBMESH_ASSERT_REALS_EQUAL( 0.5 , sum , quadrature_tolerance );
      }
  }



  //-------------------------------------------------------
  // 3D Gauss Rule Test
  template <QuadratureType qtype, Order order, unsigned int exactorder>
  void test3DWeights ()
  {
    UniquePtr<QBase> qrule = QBase::build(qtype, 3, order);
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

            LIBMESH_ASSERT_REALS_EQUAL( exact , sum , quadrature_tolerance );
          }

    // We may eventually support Gauss-Lobatto type quadrature on tets and prisms...
    if (qtype != QGAUSS_LOBATTO)
      {
        qrule->init (TET10);

        Real sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp);

        LIBMESH_ASSERT_REALS_EQUAL( 1./6., sum , quadrature_tolerance );

        qrule->init (PRISM15);

        sum = 0;

        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          sum += qrule->w(qp);

        LIBMESH_ASSERT_REALS_EQUAL( 1., sum , quadrature_tolerance );
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( QuadratureTest );
