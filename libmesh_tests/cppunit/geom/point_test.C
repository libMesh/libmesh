#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <point.h>

class PointTest : public CppUnit::TestCase { 
public: 
  CPPUNIT_TEST_SUITE( PointTest );

  CPPUNIT_TEST( testSize );
  CPPUNIT_TEST( testSizeSq );

  CPPUNIT_TEST( testValue );
  CPPUNIT_TEST( testZero );

  CPPUNIT_TEST( testEquality );
  CPPUNIT_TEST( testInEquality );
  CPPUNIT_TEST( testAssignment );

  CPPUNIT_TEST( testScalarMult );
  CPPUNIT_TEST( testScalarDiv );
  CPPUNIT_TEST( testScalarMultAssign );
  CPPUNIT_TEST( testScalarDivAssign );

  CPPUNIT_TEST( testVectorAdd );
  CPPUNIT_TEST( testVectorAddScaled );
  CPPUNIT_TEST( testVectorSub );
  CPPUNIT_TEST( testVectorMult );
  CPPUNIT_TEST( testVectorAddAssign );
  CPPUNIT_TEST( testVectorSubAssign );

  CPPUNIT_TEST_SUITE_END();

private:
  Point *m_1_1_1, *m_n1_1_n1;

public:
  void setUp()
  {
    m_1_1_1 = new Point(1,1,1);
    m_n1_1_n1 = new Point(-1,1,-1);
  }

  void tearDown() 
  {
    delete m_1_1_1;
    delete m_n1_1_n1;
  }

  void testValue()
  {
    CPPUNIT_ASSERT_EQUAL( 1.0 , (*m_1_1_1)(0));
    CPPUNIT_ASSERT_EQUAL( 1.0 , (*m_1_1_1)(1));
    CPPUNIT_ASSERT_EQUAL( 1.0 , (*m_1_1_1)(2));

    CPPUNIT_ASSERT_EQUAL( -1.0 , (*m_n1_1_n1)(0));
    CPPUNIT_ASSERT_EQUAL( 1.0  , (*m_n1_1_n1)(1));
    CPPUNIT_ASSERT_EQUAL( -1.0 , (*m_n1_1_n1)(2));
  }

  void testZero()
  {
    Point apoint(1,1,1);
    apoint.zero();
    
    CPPUNIT_ASSERT_EQUAL( 0.0 , apoint(0));
    CPPUNIT_ASSERT_EQUAL( 0.0 , apoint(1));
    CPPUNIT_ASSERT_EQUAL( 0.0 , apoint(2));
  }
  
  void testSize()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( std::sqrt(3.0) , m_1_1_1->size() , TOLERANCE*TOLERANCE );
  }

  void testSizeSq()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0 , m_1_1_1->size_sq() , TOLERANCE*TOLERANCE );
  }

  void testEquality()
  {
    CPPUNIT_ASSERT( (*m_1_1_1) == (*m_1_1_1) );    
    CPPUNIT_ASSERT( !((*m_1_1_1) == (*m_n1_1_n1)) );    
  }

  void testInEquality()
  {
    CPPUNIT_ASSERT( !((*m_1_1_1) != (*m_1_1_1)) );    
    CPPUNIT_ASSERT( (*m_1_1_1) != (*m_n1_1_n1) );    
  }

  void testAssignment()
  {
    Point apoint = (*m_1_1_1);
    
    CPPUNIT_ASSERT_EQUAL( 1.0 , (apoint)(0) );
    CPPUNIT_ASSERT_EQUAL( 1.0 , (apoint)(1) );
    CPPUNIT_ASSERT_EQUAL( 1.0 , (apoint)(2) );
  }

  void testScalarMult()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , ((*m_1_1_1)*5.0)(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , ((*m_1_1_1)*5.0)(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , ((*m_1_1_1)*5.0)(2) , TOLERANCE*TOLERANCE );
  }

  void testScalarDiv()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , ((*m_1_1_1)/5.0)(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , ((*m_1_1_1)/5.0)(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , ((*m_1_1_1)/5.0)(2) , TOLERANCE*TOLERANCE );
  }

  void testScalarMultAssign()
  {
    Point apoint(1,1,1);
    apoint*=5.0;
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , apoint(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , apoint(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , apoint(2) , TOLERANCE*TOLERANCE );
  }

  void testScalarDivAssign()
  {
    Point apoint(1.0,1.0,1.0);
    apoint/=5.0;
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , apoint(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , apoint(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , apoint(2) , TOLERANCE*TOLERANCE );
  }

  void testVectorAdd()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , ((*m_1_1_1)+(*m_n1_1_n1))(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , ((*m_1_1_1)+(*m_n1_1_n1))(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , ((*m_1_1_1)+(*m_n1_1_n1))(2) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddScaled()
  {
    Point apoint(1,1,1);
    apoint.add_scaled((*m_1_1_1),0.5);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , apoint(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , apoint(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , apoint(2) , TOLERANCE*TOLERANCE );
  }

  void testVectorSub()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , ((*m_1_1_1)-(*m_n1_1_n1))(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , ((*m_1_1_1)-(*m_n1_1_n1))(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , ((*m_1_1_1)-(*m_n1_1_n1))(2) , TOLERANCE*TOLERANCE );
  }

  void testVectorMult()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1.0 , (*m_1_1_1)*(*m_n1_1_n1) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddAssign()
  {
    Point apoint(1,1,1);
    apoint+=(*m_1_1_1);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , apoint(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , apoint(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , apoint(2) , TOLERANCE*TOLERANCE );
  }

  void testVectorSubAssign()
  {
    Point apoint(1,1,1);
    apoint-=(*m_n1_1_n1);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , apoint(0) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , apoint(1) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , apoint(2) , TOLERANCE*TOLERANCE );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( PointTest );
