// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel.h>
#include <libmesh/parallel_algebra.h>

#include "test_comm.h"

using namespace libMesh;

class ParallelPointTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelPointTest );

  CPPUNIT_TEST( testBroadcastVectorValueInt );
  CPPUNIT_TEST( testBroadcastVectorValueReal );
  CPPUNIT_TEST( testBroadcastPoint );
  CPPUNIT_TEST( testIsendRecv );
  CPPUNIT_TEST( testIrecvSend );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



  template <typename T>
  void testBroadcastVectorValue()
  {
    std::vector<VectorValue<T> > src(3), dest(3);

    {
      T val=T(0);
      for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<LIBMESH_DIM; j++)
          src[i](j) = val++;

      if (TestCommWorld->rank() == 0)
        dest = src;
    }

    TestCommWorld->broadcast(dest);

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

      if (TestCommWorld->rank() == 0)
        dest = src;
    }

    TestCommWorld->broadcast(dest);

    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<LIBMESH_DIM; j++)
        CPPUNIT_ASSERT_EQUAL (src[i](j), dest[i](j) );
  }



  void testIsendRecv ()
  {
    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    std::vector<unsigned int> src_val(3), recv_val(3);

    src_val[0] = 0;
    src_val[1] = 1;
    src_val[2] = 2;

    Parallel::Request request;

    if (TestCommWorld->size() > 1)
      {
        // Default communication
        TestCommWorld->send_mode(Parallel::Communicator::DEFAULT);

        TestCommWorld->send (procup,
                             src_val,
                             request);

        TestCommWorld->receive (procdown,
                                recv_val);

        Parallel::wait (request);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );


        // Synchronous communication
        TestCommWorld->send_mode(Parallel::Communicator::SYNCHRONOUS);
        std::fill (recv_val.begin(), recv_val.end(), 0);

        TestCommWorld->send (procup,
                             src_val,
                             request);

        TestCommWorld->receive (procdown,
                                recv_val);

        Parallel::wait (request);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );

        // Restore default communication
        TestCommWorld->send_mode(Parallel::Communicator::DEFAULT);
      }
  }



  void testIrecvSend ()
  {
    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    std::vector<unsigned int> src_val(3), recv_val(3);

    src_val[0] = 0;
    src_val[1] = 1;
    src_val[2] = 2;

    Parallel::Request request;

    if (TestCommWorld->size() > 1)
      {
        // Default communication
        TestCommWorld->send_mode(Parallel::Communicator::DEFAULT);

        TestCommWorld->receive (procdown,
                                recv_val,
                                request);

        TestCommWorld->send (procup,
                             src_val);

        Parallel::wait (request);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );

        // Synchronous communication
        TestCommWorld->send_mode(Parallel::Communicator::SYNCHRONOUS);
        std::fill (recv_val.begin(), recv_val.end(), 0);


        TestCommWorld->receive (procdown,
                                recv_val,
                                request);

        TestCommWorld->send (procup,
                             src_val);

        Parallel::wait (request);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (unsigned int i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );

        // Restore default communication
        TestCommWorld->send_mode(Parallel::Communicator::DEFAULT);
      }
  }




};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelPointTest );
