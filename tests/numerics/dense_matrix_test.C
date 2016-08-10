// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>

#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_macro.h"
#endif

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

class DenseMatrixTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(DenseMatrixTest);

  CPPUNIT_TEST(testSVD);
  CPPUNIT_TEST(testEVDreal);
  CPPUNIT_TEST(testEVDcomplex);
  CPPUNIT_TEST(testComplexSVD);

  CPPUNIT_TEST_SUITE_END();


private:

  void testSVD()
  {
    DenseMatrix<Number> U, VT;
    DenseVector<Real> sigma;
    DenseMatrix<Number> A;

    A.resize(3, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;
    A(2,0) = 5.0; A(2,1) = 6.0;

    A.svd(sigma, U, VT);

    // Solution for this case is (verified with numpy)
    DenseMatrix<Number> true_U(3,2), true_VT(2,2);
    DenseVector<Real> true_sigma(2);
    true_U(0,0) = -2.298476964000715e-01; true_U(0,1) = 8.834610176985250e-01;
    true_U(1,0) = -5.247448187602936e-01; true_U(1,1) = 2.407824921325463e-01;
    true_U(2,0) = -8.196419411205157e-01; true_U(2,1) = -4.018960334334318e-01;

    true_VT(0,0) = -6.196294838293402e-01; true_VT(0,1) = -7.848944532670524e-01;
    true_VT(1,0) = -7.848944532670524e-01; true_VT(1,1) = 6.196294838293400e-01;

    true_sigma(0) = 9.525518091565109e+00;
    true_sigma(1) = 5.143005806586446e-01;

    for (unsigned i=0; i<U.m(); ++i)
      for (unsigned j=0; j<U.n(); ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL( libmesh_real(U(i,j)), libmesh_real(true_U(i,j)), TOLERANCE*TOLERANCE);

    for (unsigned i=0; i<VT.m(); ++i)
      for (unsigned j=0; j<VT.n(); ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL( libmesh_real(VT(i,j)), libmesh_real(true_VT(i,j)), TOLERANCE*TOLERANCE);

    for (unsigned i=0; i<sigma.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma(i), true_sigma(i), TOLERANCE*TOLERANCE);
  }

  // This function is called by testEVD for different matrices.  The
  // Lapack results are compared to known eigenvalue real and
  // imaginary parts for the matrix in question, which must also be
  // passed in by non-const value, since this routine will sort them
  // in-place.
  void testEVD_helper(DenseMatrix<Real> & A,
                      std::vector<Real> true_lambda_real,
                      std::vector<Real> true_lambda_imag)
  {
    // Note: see bottom of this file, we only do this test if PETSc is
    // available, but this function currently only exists if we're
    // using real numbers.
#ifdef LIBMESH_USE_REAL_NUMBERS
    DenseVector<Real> lambda_real, lambda_imag;
    A.evd(lambda_real, lambda_imag);

    // Sort the results from Lapack *individually*.
    std::sort(lambda_real.get_values().begin(), lambda_real.get_values().end());
    std::sort(lambda_imag.get_values().begin(), lambda_imag.get_values().end());

    // Sort the true eigenvalues *individually*.
    std::sort(true_lambda_real.begin(), true_lambda_real.end());
    std::sort(true_lambda_imag.begin(), true_lambda_imag.end());

    // Compare the individually-sorted values.
    for (unsigned i=0; i<lambda_real.size(); ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(lambda_real(i), true_lambda_real[i], TOLERANCE*TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(lambda_imag(i), true_lambda_imag[i], TOLERANCE*TOLERANCE);
      }
#endif
  }

  void testEVDreal()
  {
    // This is an example from Matlab's gallery(3) which is a
    // non-symmetric 3x3 matrix with eigen values lambda = 1, 2, 3.
    DenseMatrix<Real> A(3, 3);
    A(0,0) = -149; A(0,1) = -50; A(0,2) = -154;
    A(1,0) =  537; A(1,1) = 180; A(1,2) =  546;
    A(2,0) =  -27; A(2,1) =  -9; A(2,2) =  -25;

    std::vector<Real> true_lambda_real(3);
    true_lambda_real[0] = 1.;
    true_lambda_real[1] = 2.;
    true_lambda_real[2] = 3.;
    std::vector<Real> true_lambda_imag(3); // all zero

    // call helper function to compute and verify results
    testEVD_helper(A, true_lambda_real, true_lambda_imag);
  }

  void testEVDcomplex()
  {
    // This test is also from a Matlab example, and has complex eigenvalues.
    // http://www.mathworks.com/help/matlab/math/eigenvalues.html?s_tid=gn_loc_drop
    DenseMatrix<Real> A(3, 3);
    A(0,0) =  0; A(0,1) = -6; A(0,2) =  -1;
    A(1,0) =  6; A(1,1) =  2; A(1,2) = -16;
    A(2,0) = -5; A(2,1) = 20; A(2,2) = -10;

    std::vector<Real> true_lambda_real(3);
    true_lambda_real[0] = -3.070950351248293;
    true_lambda_real[1] = -2.464524824375853;
    true_lambda_real[2] = -2.464524824375853;
    std::vector<Real> true_lambda_imag(3);
    true_lambda_imag[0] = 0.;
    true_lambda_imag[1] = 17.60083096447099;
    true_lambda_imag[2] = -17.60083096447099;

    // call helper function to compute and verify results
    testEVD_helper(A, true_lambda_real, true_lambda_imag);
  }

  void testComplexSVD()
  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    DenseMatrix<Complex> A(3,3);

    A(0,0) = Complex(2.18904,4.44523e-18); A(0,1) = Complex(-3.20491,-0.136699);   A(0,2) = Complex(0.716316,-0.964802);
    A(1,0) = Complex(-3.20491,0.136699);   A(1,1) = Complex(4.70076,-3.25261e-18); A(1,2) = Complex(-0.98849,1.45727);
    A(2,0) = Complex(0.716316,0.964802);   A(2,1) = Complex(-0.98849,-1.45727);    A(2,2) = Complex(0.659629,-4.01155e-18);

    DenseVector<Real> sigma;
    A.svd(sigma);

    DenseVector<Real> true_sigma(3);
    true_sigma(0) = 7.54942516052;
    true_sigma(1) = 3.17479511368e-06;
    true_sigma(2) = 6.64680908281e-07;

    for (unsigned i=0; i<sigma.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma(i), true_sigma(i), 1.e-10);
#endif
  }

};

// Only run the test if we expect it can actually work!
#ifdef LIBMESH_HAVE_PETSC
#if !PETSC_VERSION_LESS_THAN(3,1,0)
  CPPUNIT_TEST_SUITE_REGISTRATION(DenseMatrixTest);
#endif
#endif
