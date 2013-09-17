#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <libmesh/parallel.h>
#include <libmesh/parallel_algebra.h>

using namespace libMesh;

class ParallelTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelTest );

#ifndef LIBMESH_DISABLE_COMMWORLD
  CPPUNIT_TEST( testGather );
  CPPUNIT_TEST( testAllGather );
  CPPUNIT_TEST( testBroadcast );
  CPPUNIT_TEST( testBroadcastVectorValueInt );
  CPPUNIT_TEST( testBroadcastVectorValueReal );
  CPPUNIT_TEST( testBroadcastPoint );
  CPPUNIT_TEST( testBarrier );
  CPPUNIT_TEST( testMin );
  CPPUNIT_TEST( testMax );
  CPPUNIT_TEST( testIsendRecv );
  CPPUNIT_TEST( testIrecvSend );
#endif // !LIBMESH_DISABLE_COMMWORLD

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



#ifndef LIBMESH_DISABLE_COMMWORLD
  void testGather()
  {
    std::vector<processor_id_type> vals;
    CommWorld.gather(0,libMesh::processor_id(),vals);

    if (libMesh::processor_id() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
	CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



  void testAllGather()
  {
    std::vector<processor_id_type> vals;
    CommWorld.allgather(libMesh::processor_id(),vals);

    for (processor_id_type i=0; i<vals.size(); i++)
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

    CommWorld.broadcast(dest);

    for (unsigned int i=0; i<src.size(); i++)
      CPPUNIT_ASSERT_EQUAL( src[i] , dest[i] );
  }



  template <typename T>
  void testBroadcastVectorValue()
  {
    std::vector<VectorValue<T> > src(3), dest(3);

    {
      T val=T(0);
      for (unsigned int i=0; i<3; i++)
	for (unsigned int j=0; j<LIBMESH_DIM; j++)
	  src[i](j) = val++;

      if (libMesh::processor_id() == 0)
	dest = src;
    }

    CommWorld.broadcast(dest);

    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<LIBMESH_DIM; j++)
	CPPUNIT_ASSERT_EQUAL (src[i](j), dest[i](j) );
  }



  void testBroadcastVectorValueInt()
  {
    this->testBroadcastVectorValue<int>();
  }



  void testBroadcastVectorValueReal()
  {
    this->testBroadcastVectorValue<Real>();
  }



  void testBroadcastPoint()
  {
    std::vector<Point> src(3), dest(3);

    {
      Real val=0.;
      for (unsigned int i=0; i<3; i++)
	for (unsigned int j=0; j<LIBMESH_DIM; j++)
	  src[i](j) = val++;

      if (libMesh::processor_id() == 0)
	dest = src;
    }

    CommWorld.broadcast(dest);

    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<LIBMESH_DIM; j++)
	CPPUNIT_ASSERT_EQUAL (src[i](j), dest[i](j) );
  }



  void testBarrier()
  {
    CommWorld.barrier();
  }



  void testMin ()
  {
    unsigned int min = libMesh::processor_id();

    CommWorld.min(min);

    CPPUNIT_ASSERT_EQUAL (min, static_cast<unsigned int>(0));
  }



  void testMax ()
  {
    processor_id_type max = libMesh::processor_id();

    CommWorld.max(max);

    CPPUNIT_ASSERT_EQUAL (static_cast<processor_id_type>(max+1),
			  libMesh::n_processors());
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
        CommWorld.send (procup,
		        src_val,
		        request);

        CommWorld.receive (procdown,
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
        CommWorld.receive (procdown,
		           recv_val,
		           request);

        CommWorld.send (procup,
		        src_val);

        Parallel::wait (request);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );
      }
  }
#endif // !LIBMESH_DISABLE_COMMWORLD

};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelTest );
