#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <quadrature.h>

class QuadratureTest : public CppUnit::TestCase { 
public: 
  CPPUNIT_TEST_SUITE( QuadratureTest );
  
  CPPUNIT_TEST( testBuild );
  CPPUNIT_TEST( testWeights );

  CPPUNIT_TEST_SUITE_END();

private:


public:
  void setUp ()
  {}

  void tearDown ()
  {}



  void testBuild ()
  {
    AutoPtr<QBase> qrule1D = QBase::build (QGAUSS, 1, THIRD);
    AutoPtr<QBase> qrule2D = QBase::build (QGAUSS, 2, THIRD);
    AutoPtr<QBase> qrule3D = QBase::build (QGAUSS, 3, THIRD);

    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(1) , qrule1D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(2) , qrule2D->get_dim() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned int>(3) , qrule3D->get_dim() );
  }



  void testWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, 3, FIFTH);
    qrule->init (HEX8);

    Real sum = 0;

    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0 , sum , TOLERANCE*TOLERANCE );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( QuadratureTest );
