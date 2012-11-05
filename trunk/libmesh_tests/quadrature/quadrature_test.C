#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <libmesh/quadrature.h>

using namespace libMesh;

class QuadratureTest : public CppUnit::TestCase { 
public: 
  CPPUNIT_TEST_SUITE( QuadratureTest );
  
  CPPUNIT_TEST( testBuild );
  
  CPPUNIT_TEST( test1DWeights<FIRST> );
  CPPUNIT_TEST( test1DWeights<SECOND> );
  CPPUNIT_TEST( test1DWeights<THIRD> );
  CPPUNIT_TEST( test1DWeights<FOURTH> );
  CPPUNIT_TEST( test1DWeights<FIFTH> );
  CPPUNIT_TEST( test1DWeights<SIXTH> );
  CPPUNIT_TEST( test1DWeights<SEVENTH> );

  CPPUNIT_TEST( test2DWeights<FIRST> );
  CPPUNIT_TEST( test2DWeights<SECOND> );
  CPPUNIT_TEST( test2DWeights<THIRD> );
  CPPUNIT_TEST( test2DWeights<FOURTH> );
  CPPUNIT_TEST( test2DWeights<FIFTH> );
  CPPUNIT_TEST( test2DWeights<SIXTH> );
  CPPUNIT_TEST( test2DWeights<SEVENTH> );

  CPPUNIT_TEST( test3DWeights<FIRST> );
  CPPUNIT_TEST( test3DWeights<SECOND> );
  CPPUNIT_TEST( test3DWeights<THIRD> );
  CPPUNIT_TEST( test3DWeights<FOURTH> );
  CPPUNIT_TEST( test3DWeights<FIFTH> );
  CPPUNIT_TEST( test3DWeights<SIXTH> );
  CPPUNIT_TEST( test3DWeights<SEVENTH> );

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



  //-------------------------------------------------------
  // 1D Gauss Rule Test
  template <Order order>
  void test1DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, 1, order);
    qrule->init (EDGE3);
    
    Real sum = 0;
    
    for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      sum += qrule->w(qp);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , sum , TOLERANCE*TOLERANCE );
  }
    


  //-------------------------------------------------------
  // 2D Gauss Rule Test
  template <Order order>
  void test2DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, 2, order);
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
  template <Order order>
  void test3DWeights ()
  {
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, 3, order);
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
