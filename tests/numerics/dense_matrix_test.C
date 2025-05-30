// libmesh includes
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>

#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_macro.h"
#endif

#include "libmesh_cppunit.h"


using namespace libMesh;

class DenseMatrixTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  LIBMESH_CPPUNIT_TEST_SUITE(DenseMatrixTest);

  CPPUNIT_TEST(testOuterProduct);
  CPPUNIT_TEST(testSVD);
  CPPUNIT_TEST(testEVDreal);
  CPPUNIT_TEST(testEVDcomplex);
  CPPUNIT_TEST(testComplexSVD);
  CPPUNIT_TEST(testSubMatrix);

  CPPUNIT_TEST_SUITE_END();


private:

  void testOuterProduct()
  {
    LOG_UNIT_TEST;

    DenseVector<Real> a = {1.0, 2.0};
    DenseVector<Real> b = {3.0, 4.0, 5.0};

    DenseMatrix<Real> a_times_b;
    a_times_b.outer_product(a, b);

    DenseMatrix<Real> a_times_b_correct(2, 3,
                                        {3., 4., 5.,
                                         6., 8., 10.});

    for (unsigned int i = 0; i < a.size(); ++i)
      for (unsigned int j = 0; j < b.size(); ++j)
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (a_times_b(i,j), a_times_b_correct(i,j), TOLERANCE*TOLERANCE);
  }

  void testSVD()
  {
    LOG_UNIT_TEST;

    DenseMatrix<Number> U, VT;
    DenseVector<Real> sigma;
    DenseMatrix<Number> A(3, 2, {
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0});

    A.svd(sigma, U, VT);

    // Solution for this case is (verified with numpy)
    DenseMatrix<Number> true_U(3, 2, {
        -2.298476964000715e-01, 8.834610176985250e-01,
        -5.247448187602936e-01, 2.407824921325463e-01,
        -8.196419411205157e-01, -4.018960334334318e-01});

    DenseMatrix<Number> true_VT(2, 2, {
        -6.196294838293402e-01, -7.848944532670524e-01,
        -7.848944532670524e-01, 6.196294838293400e-01});

    DenseVector<Real> true_sigma = {9.525518091565109e+00, 5.143005806586446e-01};

    // Tolerance is bounded by the double literals above and by using Real
    const Real tol = std::max(Real(1e-12), TOLERANCE*TOLERANCE);

    for (unsigned i=0; i<U.m(); ++i)
      for (unsigned j=0; j<U.n(); ++j)
        LIBMESH_ASSERT_NUMBERS_EQUAL( U(i,j), true_U(i,j), tol);

    for (unsigned i=0; i<VT.m(); ++i)
      for (unsigned j=0; j<VT.n(); ++j)
        LIBMESH_ASSERT_NUMBERS_EQUAL( VT(i,j), true_VT(i,j), tol);

    for (unsigned i=0; i<sigma.size(); ++i)
      LIBMESH_ASSERT_FP_EQUAL(sigma(i), true_sigma(i), tol);
  }

  // This function is called by testEVD for different matrices.  The
  // Lapack results are compared to known eigenvalue real and
  // imaginary parts for the matrix in question, which must also be
  // passed in by non-const value, since this routine will sort them
  // in-place.
  void testEVD_helper(DenseMatrix<Real> & A,
                      std::vector<Real> true_lambda_real,
                      std::vector<Real> true_lambda_imag,
                      Real tol)
  {
    // Note: see bottom of this file, we only do this test if PETSc is
    // available, but this function currently only exists if we're
    // using real numbers.
#ifdef LIBMESH_USE_REAL_NUMBERS
    // Let's compute the eigenvalues on a copy of A, so that we can
    // use the original to check the computation.
    DenseMatrix<Real> A_copy = A;

    DenseVector<Real> lambda_real, lambda_imag;
    DenseMatrix<Real> VR; // right eigenvectors
    DenseMatrix<Real> VL; // left eigenvectors
    A_copy.evd_left_and_right(lambda_real, lambda_imag, VL, VR);

    // The matrix is square and of size N x N.
    const unsigned N = A.m();

    // Verify left eigen-values.
    // Test that the right eigenvalues are self-consistent by computing
    // u_j**H * A = lambda_j * u_j**H
    // Note that we have to handle real and complex eigenvalues
    // differently, since complex eigenvectors share their storage.
    for (unsigned eigenval=0; eigenval<N; ++eigenval)
      {
        // Only check real eigenvalues
        if (std::abs(lambda_imag(eigenval)) < tol*tol)
          {
            // remove print libMesh::out << "Checking eigenvalue: " << eigenval << std::endl;
            DenseVector<Real> lhs(N), rhs(N);
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = lambda_real(eigenval) * VL(i, eigenval);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(j, i) * VL(j, eigenval); // Note: A(j,i)
              }

            // Subtract and assert that the norm of the difference is
            // below some tolerance.
            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);
          }
        else
          {
            // This is a complex eigenvalue, so:
            // a.) It occurs in a complex-conjugate pair
            // b.) the real part of the eigenvector is stored is VL(:,eigenval)
            // c.) the imag part of the eigenvector is stored in VL(:,eigenval+1)
            //
            // Equating the real and imaginary parts of Ax=lambda*x leads to two sets
            // of relations that must hold:
            // 1.) A^T x_r =  lambda_r*x_r + lambda_i*x_i
            // 2.) A^T x_i = -lambda_i*x_r + lambda_r*x_i
            // which we can verify.

            // 1.)
            DenseVector<Real> lhs(N), rhs(N);
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = lambda_real(eigenval) * VL(i, eigenval) + lambda_imag(eigenval) * VL(i, eigenval+1);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(j, i) * VL(j, eigenval); // Note: A(j,i)
              }

            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);

            // libMesh::out << "lhs=" << std::endl;
            // lhs.print_scientific(libMesh::out, /*precision=*/15);
            //
            // libMesh::out << "rhs=" << std::endl;
            // rhs.print_scientific(libMesh::out, /*precision=*/15);

            // 2.)
            lhs.zero();
            rhs.zero();
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = -lambda_imag(eigenval) * VL(i, eigenval) + lambda_real(eigenval) * VL(i, eigenval+1);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(j, i) * VL(j, eigenval+1); // Note: A(j,i)
              }

            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);

            // libMesh::out << "lhs=" << std::endl;
            // lhs.print_scientific(libMesh::out, /*precision=*/15);
            //
            // libMesh::out << "rhs=" << std::endl;
            // rhs.print_scientific(libMesh::out, /*precision=*/15);

            // We'll skip the second member of the complex conjugate
            // pair.  If the first one worked, the second one should
            // as well...
            eigenval += 1;
          }
      }

    // Verify right eigen-values.
    // Test that the right eigenvalues are self-consistent by computing
    // A * v_j - lambda_j * v_j
    // Note that we have to handle real and complex eigenvalues
    // differently, since complex eigenvectors share their storage.
    for (unsigned eigenval=0; eigenval<N; ++eigenval)
      {
        // Only check real eigenvalues
        if (std::abs(lambda_imag(eigenval)) < tol*tol)
          {
            // remove print libMesh::out << "Checking eigenvalue: " << eigenval << std::endl;
            DenseVector<Real> lhs(N), rhs(N);
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = lambda_real(eigenval) * VR(i, eigenval);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(i, j) * VR(j, eigenval);
              }

            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);
          }
        else
          {
            // This is a complex eigenvalue, so:
            // a.) It occurs in a complex-conjugate pair
            // b.) the real part of the eigenvector is stored is VR(:,eigenval)
            // c.) the imag part of the eigenvector is stored in VR(:,eigenval+1)
            //
            // Equating the real and imaginary parts of Ax=lambda*x leads to two sets
            // of relations that must hold:
            // 1.) Ax_r = lambda_r*x_r - lambda_i*x_i
            // 2.) Ax_i = lambda_i*x_r + lambda_r*x_i
            // which we can verify.

            // 1.)
            DenseVector<Real> lhs(N), rhs(N);
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = lambda_real(eigenval) * VR(i, eigenval) - lambda_imag(eigenval) * VR(i, eigenval+1);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(i, j) * VR(j, eigenval);
              }

            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);

            // 2.)
            lhs.zero();
            rhs.zero();
            for (unsigned i=0; i<N; ++i)
              {
                rhs(i) = lambda_imag(eigenval) * VR(i, eigenval) + lambda_real(eigenval) * VR(i, eigenval+1);
                for (unsigned j=0; j<N; ++j)
                  lhs(i) += A(i, j) * VR(j, eigenval+1);
              }

            lhs -= rhs;
            LIBMESH_ASSERT_FP_EQUAL(/*expected=*/0., /*actual=*/lhs.l2_norm(), std::sqrt(tol)*tol);

            // We'll skip the second member of the complex conjugate
            // pair.  If the first one worked, the second one should
            // as well...
            eigenval += 1;
          }
      }

    // Sort the results from Lapack *individually*.
    std::sort(lambda_real.get_values().begin(), lambda_real.get_values().end());
    std::sort(lambda_imag.get_values().begin(), lambda_imag.get_values().end());

    // Sort the true eigenvalues *individually*.
    std::sort(true_lambda_real.begin(), true_lambda_real.end());
    std::sort(true_lambda_imag.begin(), true_lambda_imag.end());

    // Compare the individually-sorted values.
    for (unsigned i=0; i<lambda_real.size(); ++i)
      {
        // Note: I initially verified the results with TOLERANCE**2,
        // but that turned out to be just a bit too tight for some of
        // the test problems.  I'm not sure what controls the accuracy
        // of the eigenvalue computation in LAPACK, there is no way to
        // set a tolerance in the LAPACKgeev_ interface.
        LIBMESH_ASSERT_FP_EQUAL(/*expected=*/true_lambda_real[i], /*actual=*/lambda_real(i), std::sqrt(tol)*tol);
        LIBMESH_ASSERT_FP_EQUAL(/*expected=*/true_lambda_imag[i], /*actual=*/lambda_imag(i), std::sqrt(tol)*tol);
      }
#endif
  }

  void testEVDreal()
  {
    LOG_UNIT_TEST;

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

    // call helper function to compute and verify results.
    testEVD_helper(A, true_lambda_real, true_lambda_imag, TOLERANCE);
  }

  void testEVDcomplex()
  {
    LOG_UNIT_TEST;

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

    // We're good up to double precision (due to the truncated
    // literals above) or Real precision, whichever is worse
    const Real tol = std::max(Real(1e-6), TOLERANCE);

    // call helper function to compute and verify results
    testEVD_helper(A, true_lambda_real, true_lambda_imag, tol);
  }

  void testComplexSVD()
  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    LOG_UNIT_TEST;

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
      LIBMESH_ASSERT_FP_EQUAL(sigma(i), true_sigma(i), 1.e-10);
#endif
  }

  void testSubMatrix()
  {
    LOG_UNIT_TEST;

    DenseMatrix<Number> A(4, 3);
    A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;
    A(2,0) = 7.0; A(2,1) = 8.0; A(2,2) = 9.0;
    A(3,0) =10.0; A(3,1) =11.0; A(3,2) =12.0;

    DenseMatrix<Number> B(2, 2);
    B(0,0) = 7.0; B(0,1) = 8.0;
    B(1,0) =10.0; B(1,1) =11.0;

    DenseMatrix<Number> C = A.sub_matrix(2, 2, 0, 2);
    CPPUNIT_ASSERT(B == C);
  }
};

// These tests require PETSc
#ifdef LIBMESH_HAVE_PETSC
CPPUNIT_TEST_SUITE_REGISTRATION(DenseMatrixTest);
#endif
