#ifndef __type_vector_test_h__
#define __type_vector_test_h__

#include <libmesh/type_vector.h>

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#define TYPEVECTORTEST                          \
  CPPUNIT_TEST( testNorm );                     \
  CPPUNIT_TEST( testNormSq );                   \
                                                \
  CPPUNIT_TEST( testValue );                    \
  CPPUNIT_TEST( testZero );                     \
                                                \
  CPPUNIT_TEST( testEquality );                 \
  CPPUNIT_TEST( testInEquality );               \
  CPPUNIT_TEST( testAssignment );               \
                                                \
  CPPUNIT_TEST( testScalarMult );               \
  CPPUNIT_TEST( testScalarDiv );                \
  CPPUNIT_TEST( testScalarMultAssign );         \
  CPPUNIT_TEST( testScalarDivAssign );          \
                                                \
  CPPUNIT_TEST( testVectorAdd );                \
  CPPUNIT_TEST( testVectorAddScaled );          \
  CPPUNIT_TEST( testVectorSub );                \
  CPPUNIT_TEST( testVectorMult );               \
  CPPUNIT_TEST( testVectorAddAssign );          \
  CPPUNIT_TEST( testVectorSubAssign );          \
                                                \
  CPPUNIT_TEST( testNormBase );                 \
  CPPUNIT_TEST( testNormSqBase );               \
  CPPUNIT_TEST( testValueBase );                \
  CPPUNIT_TEST( testZeroBase );                 \
                                                \
  CPPUNIT_TEST( testEqualityBase );             \
  CPPUNIT_TEST( testInEqualityBase );           \
  CPPUNIT_TEST( testAssignmentBase );           \
                                                \
  CPPUNIT_TEST( testScalarMultBase );           \
  CPPUNIT_TEST( testScalarDivBase );            \
  CPPUNIT_TEST( testScalarMultAssignBase );     \
  CPPUNIT_TEST( testScalarDivAssignBase );      \
                                                \
  CPPUNIT_TEST( testVectorAddBase );            \
  CPPUNIT_TEST( testVectorAddScaledBase );      \
  CPPUNIT_TEST( testVectorSubBase );            \
  CPPUNIT_TEST( testVectorMultBase );           \
  CPPUNIT_TEST( testVectorAddAssignBase );      \
  CPPUNIT_TEST( testVectorSubAssignBase );


using namespace libMesh;

template <class DerivedClass>
class TypeVectorTestBase : public CppUnit::TestCase {

protected:
  typedef typename DerivedClass::value_type T;

private:
  DerivedClass *m_1_1_1, *m_n1_1_n1;
  TypeVector<T>   *basem_1_1_1, *basem_n1_1_n1;

public:
  virtual void setUp()
  {
    m_1_1_1 = new DerivedClass(1,1,1);
    m_n1_1_n1 = new DerivedClass(-1,1,-1);

    basem_1_1_1 = m_1_1_1;
    basem_n1_1_n1 = m_n1_1_n1;
  }

  virtual void tearDown()
  {
    delete m_1_1_1;
    delete m_n1_1_n1;
  }

  void testValue()
  {
    CPPUNIT_ASSERT_EQUAL( T(1), (*m_1_1_1)(0));
    CPPUNIT_ASSERT_EQUAL( T(1), (*m_1_1_1)(1));
    CPPUNIT_ASSERT_EQUAL( T(1), (*m_1_1_1)(2));

    CPPUNIT_ASSERT_EQUAL( T(-1), (*m_n1_1_n1)(0));
    CPPUNIT_ASSERT_EQUAL( T(1) , (*m_n1_1_n1)(1));
    CPPUNIT_ASSERT_EQUAL( T(-1), (*m_n1_1_n1)(2));
  }

  void testZero()
  {
    DerivedClass avector(1,1,1);
    avector.zero();

    CPPUNIT_ASSERT_EQUAL( T(0), avector(0));
    CPPUNIT_ASSERT_EQUAL( T(0), avector(1));
    CPPUNIT_ASSERT_EQUAL( T(0), avector(2));
  }

  void testNorm()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( std::sqrt(3.0) , m_1_1_1->norm() , TOLERANCE*TOLERANCE );
  }

  void testNormSq()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0 , m_1_1_1->norm_sq() , TOLERANCE*TOLERANCE );
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
    DerivedClass avector = (*m_1_1_1);

    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(0) );
    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(1) );
    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(2) );
  }

  void testScalarInit()  // Not valid for all derived classes!
  {
    DerivedClass avector = 0;
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );

    DerivedClass bvector = 2.0;
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(bvector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(bvector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(bvector(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarMult()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*m_1_1_1)*5.0)(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*m_1_1_1)*5.0)(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*m_1_1_1)*5.0)(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarDiv()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*m_1_1_1)/5.0)(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*m_1_1_1)/5.0)(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*m_1_1_1)/5.0)(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarMultAssign()
  {
    DerivedClass avector(1,1,1);
    avector*=5.0;

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarDivAssign()
  {
    DerivedClass avector(1.0,1.0,1.0);
    avector/=5.0;

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAdd()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*m_1_1_1)+(*m_n1_1_n1))(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*m_1_1_1)+(*m_n1_1_n1))(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*m_1_1_1)+(*m_n1_1_n1))(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddScaled()
  {
    DerivedClass avector(1,1,1);
    avector.add_scaled((*m_1_1_1),0.5);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorSub()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*m_1_1_1)-(*m_n1_1_n1))(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*m_1_1_1)-(*m_n1_1_n1))(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*m_1_1_1)-(*m_n1_1_n1))(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorMult()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1.0 , libmesh_real((*m_1_1_1)*(*m_n1_1_n1)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddAssign()
  {
    DerivedClass avector(1,1,1);
    avector+=(*m_1_1_1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorSubAssign()
  {
    DerivedClass avector(1,1,1);
    avector-=(*m_n1_1_n1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testValueBase()
  {
    CPPUNIT_ASSERT_EQUAL( T(1), (*basem_1_1_1)(0));
    CPPUNIT_ASSERT_EQUAL( T(1), (*basem_1_1_1)(1));
    CPPUNIT_ASSERT_EQUAL( T(1), (*basem_1_1_1)(2));

    CPPUNIT_ASSERT_EQUAL( T(-1), (*basem_n1_1_n1)(0));
    CPPUNIT_ASSERT_EQUAL( T(1 ), (*basem_n1_1_n1)(1));
    CPPUNIT_ASSERT_EQUAL( T(-1), (*basem_n1_1_n1)(2));
  }

  void testZeroBase()
  {
    TypeVector<T> avector((*basem_1_1_1));
    avector.zero();

    CPPUNIT_ASSERT_EQUAL( T(0), avector(0));
    CPPUNIT_ASSERT_EQUAL( T(0), avector(1));
    CPPUNIT_ASSERT_EQUAL( T(0), avector(2));
  }

  void testNormBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( std::sqrt(3.0) , basem_1_1_1->norm() , TOLERANCE*TOLERANCE );
  }

  void testNormSqBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0 , basem_1_1_1->norm_sq() , TOLERANCE*TOLERANCE );
  }

  void testEqualityBase()
  {
    CPPUNIT_ASSERT( (*basem_1_1_1) == (*basem_1_1_1) );
    CPPUNIT_ASSERT( !((*basem_1_1_1) == (*basem_n1_1_n1)) );
  }

  void testInEqualityBase()
  {
    CPPUNIT_ASSERT( !((*basem_1_1_1) != (*basem_1_1_1)) );
    CPPUNIT_ASSERT( (*basem_1_1_1) != (*basem_n1_1_n1) );
  }

  void testAssignmentBase()
  {
    TypeVector<T> avector = (*m_1_1_1);

    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(0) );
    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(1) );
    CPPUNIT_ASSERT_EQUAL( T(1), (avector)(2) );
  }

  void testScalarMultBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*basem_1_1_1)*5.0)(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*basem_1_1_1)*5.0)(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(((*basem_1_1_1)*5.0)(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarDivBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*basem_1_1_1)/5.0)(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*basem_1_1_1)/5.0)(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(((*basem_1_1_1)/5.0)(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarMultAssignBase()
  {
    TypeVector<T> avector(*m_1_1_1);
    avector*=5.0;

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testScalarDivAssignBase()
  {
    TypeVector<T> avector(*m_1_1_1);
    avector/=5.0;

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0/5.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*basem_1_1_1)+(*basem_n1_1_n1))(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*basem_1_1_1)+(*basem_n1_1_n1))(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*basem_1_1_1)+(*basem_n1_1_n1))(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddScaledBase()
  {
    TypeVector<T> avector(*m_1_1_1);
    avector.add_scaled((*basem_1_1_1),0.5);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorSubBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*basem_1_1_1)-(*basem_n1_1_n1))(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(((*basem_1_1_1)-(*basem_n1_1_n1))(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(((*basem_1_1_1)-(*basem_n1_1_n1))(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorMultBase()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1.0 , libmesh_real((*basem_1_1_1)*(*basem_n1_1_n1)) , TOLERANCE*TOLERANCE );
  }

  void testVectorAddAssignBase()
  {
    TypeVector<T> avector(*m_1_1_1);
    avector+=(*basem_1_1_1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }

  void testVectorSubAssignBase()
  {
    TypeVector<T> avector(*m_1_1_1);
    avector-=(*basem_n1_1_n1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(0)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0 , libmesh_real(avector(1)) , TOLERANCE*TOLERANCE );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0 , libmesh_real(avector(2)) , TOLERANCE*TOLERANCE );
  }
};

#endif // #ifdef __type_vector_test_h__
