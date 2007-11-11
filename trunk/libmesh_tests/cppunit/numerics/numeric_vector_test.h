#ifndef __numeric_vector_test_h__
#define __numeric_vector_test_h__

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#define NUMERICVECTORTEST \
  CPPUNIT_TEST( testLocalize ); \
  CPPUNIT_TEST( testLocalizeBase ); \
  CPPUNIT_TEST( testLocalizeToOne ); \
  CPPUNIT_TEST( testLocalizeToOneBase );

template <class DerivedClass>
class NumericVectorTest : public CppUnit::TestCase { 

public:
  void setUp()
  {}

  void tearDown() 
  {}

  template <class Base, class Derived>
  void Localize(bool to_one=false)
  {
    const unsigned int root_pid = 0;
    unsigned int block_size  = 10;
    unsigned int local_size  = block_size + libMesh::processor_id();  // a different size on
    unsigned int global_size = 0;                                     // each processor.
    
    for (unsigned int p=0; p<libMesh::n_processors(); p++)
      global_size += (block_size + p);
    
    {
      Base & v = *(new Derived(global_size, local_size));
      std::vector<Real> l(global_size);
      
      const unsigned int 
	first = v.first_local_index(),
	last  = v.last_local_index();
      
      for (unsigned int n=first; n != last; n++)
	v.set (n, static_cast<Real>(n));
      v.close();
      
      if(!to_one)
	v.localize(l);
      else
	v.localize_to_one(l,root_pid);

      if(!to_one || libMesh::processor_id() == root_pid)
	//Yes I really mean v.size()
	for (unsigned int i=0; i<v.size(); i++)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL( i , l[i] , TOLERANCE*TOLERANCE );

      delete &v;
    }
  }
  
  void testLocalize()
  {
    Localize<DerivedClass,DerivedClass>();
  }

  void testLocalizeBase()
  {
    Localize<NumericVector<Real>,DerivedClass>();
  }

  void testLocalizeToOne()
  {
    Localize<DerivedClass,DerivedClass >(true);
  }

  void testLocalizeToOneBase()
  {
    Localize<NumericVector<Real>,DerivedClass>(true);
  }
};

#endif // #ifdef __numeric_vector_test_h__
