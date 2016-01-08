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

#ifdef LIBMESH_USE_REAL_NUMBERS

using namespace libMesh;

class DenseMatrixTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(DenseMatrixTest);

  CPPUNIT_TEST(testSVD);

  CPPUNIT_TEST_SUITE_END();


private:
  void testSVD()
  {
    DenseMatrix<Number> U, VT;
    DenseVector<Number> sigma;
    DenseMatrix<Number> A;

    A.resize(3, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;
    A(2,0) = 5.0; A(2,1) = 6.0;

    A.svd(sigma, U, VT);

    // Solution for this case is (verified with numpy)
    DenseMatrix<Number> true_U(3,2), true_VT(2,2);
    DenseVector<Number> true_sigma(2);
    true_U(0,0) = -2.298476964000715e-01; true_U(0,1) = 8.834610176985250e-01;
    true_U(1,0) = -5.247448187602936e-01; true_U(1,1) = 2.407824921325463e-01;
    true_U(2,0) = -8.196419411205157e-01; true_U(2,1) = -4.018960334334318e-01;

    true_VT(0,0) = -6.196294838293402e-01; true_VT(0,1) = -7.848944532670524e-01;
    true_VT(1,0) = -7.848944532670524e-01; true_VT(1,1) = 6.196294838293400e-01;

    true_sigma(0) = 9.525518091565109e+00;
    true_sigma(1) = 5.143005806586446e-01;

    for (unsigned i=0; i<U.m(); ++i)
      for (unsigned j=0; j<U.n(); ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(U(i,j), true_U(i,j), 1.e-12);

    for (unsigned i=0; i<VT.m(); ++i)
      for (unsigned j=0; j<VT.n(); ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(VT(i,j), true_VT(i,j), 1.e-12);

    for (unsigned i=0; i<sigma.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma(i), true_sigma(i), 1.e-12);
  }

};

// Only run the test if we expect it can actually work!
#ifdef LIBMESH_HAVE_PETSC
#if !PETSC_VERSION_LESS_THAN(3,1,0)
  CPPUNIT_TEST_SUITE_REGISTRATION(DenseMatrixTest);
#endif
#endif

#endif // LIBMESH_USE_REAL_NUMBERS
