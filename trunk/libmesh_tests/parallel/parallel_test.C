#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <libmesh/parallel.h>

class ParallelTest : public CppUnit::TestCase { 
public: 
  CPPUNIT_TEST_SUITE( ParallelTest );

  CPPUNIT_TEST( testGather );
  CPPUNIT_TEST( testAllGather );
  CPPUNIT_TEST( testBroadcast );
  CPPUNIT_TEST( testBarrier );
  CPPUNIT_TEST( testMin );
  CPPUNIT_TEST( testMax );
  CPPUNIT_TEST( testIsendRecv );
  CPPUNIT_TEST( testIrecvSend );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown() 
  {}



  void testGather()
  {
    std::vector<unsigned int> vals;
    Parallel::gather(0,libMesh::processor_id(),vals);
    
    if (libMesh::processor_id() == 0)
      for (unsigned int i=0; i<vals.size(); i++)
	CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



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



  void testBarrier()
  {
    Parallel::barrier();
  }



  void testMin ()
  {
    unsigned int min = libMesh::processor_id();
    
    Parallel::min(min);
    
    CPPUNIT_ASSERT_EQUAL (min, static_cast<unsigned int>(0));
  }



  void testMax ()
  {
    unsigned int max = libMesh::processor_id();
    
    Parallel::max(max);
    
    CPPUNIT_ASSERT_EQUAL (max+1, libMesh::n_processors());
  }



  void testIsendRecv ()
  {
    unsigned int procup = (libMesh::processor_id() + 1) %
                           libMesh::n_processors();
    unsigned int procdown = (libMesh::n_processors() +
			     libMesh::processor_id() - 1) %
                             libMesh::n_processors();
    
    std::vector<unsigned int> src_val(3), recv_val(3);

    src_val[0] = 0;
    src_val[1] = 1;
    src_val[2] = 2;

    Parallel::Request request;

    if (libMesh::n_processors() > 1)
      {
        Parallel::send (procup,
		        src_val,
		        request);

        Parallel::receive (procdown,
		           recv_val);

        Parallel::wait (request);
    
        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );
      }
  }



  void testIrecvSend ()
  {
    unsigned int procup = (libMesh::processor_id() + 1) %
                           libMesh::n_processors();
    unsigned int procdown = (libMesh::n_processors() +
			     libMesh::processor_id() - 1) %
                             libMesh::n_processors();
    
    std::vector<unsigned int> src_val(3), recv_val(3);

    src_val[0] = 0;
    src_val[1] = 1;
    src_val[2] = 2;

    Parallel::Request request;

    if (libMesh::n_processors() > 1)
      {
        Parallel::receive (procdown,
		           recv_val,
		           request);

        Parallel::send (procup,
		        src_val);

        Parallel::wait (request);
    
        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelTest );
