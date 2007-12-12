#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <parallel.h>

class ParallelTest : public CppUnit::TestCase { 
public: 
  CPPUNIT_TEST_SUITE( ParallelTest );

  CPPUNIT_TEST( testAllGather );
  CPPUNIT_TEST( testBroadcast );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown() 
  {}

  void testAllGather()
  {
    std::vector<unsigned int> vals;
    Parallel::allgather(libMesh::processor_id(),vals);
    
    for (unsigned int i=0; i<vals.size(); i++)
      CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }

  void testBroadcast()
  {
    std::vector<unsigned int> src(3), dest(3);

    src[0]=0;
    src[1]=1;
    src[2]=2;

    if (libMesh::processor_id() == 0)
      dest = src;

    Parallel::broadcast(dest);

    for (unsigned int i=0; i<src.size(); i++)
      CPPUNIT_ASSERT_EQUAL( src[i] , dest[i] );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelTest );
