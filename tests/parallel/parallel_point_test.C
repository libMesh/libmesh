// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel.h>
#include <libmesh/parallel_algebra.h>

#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

class ParallelPointTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelPointTest );

  CPPUNIT_TEST( testAllGatherPoint );
  CPPUNIT_TEST( testAllGatherPairPointPoint );
  CPPUNIT_TEST( testAllGatherPairRealPoint );
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



  void testAllGatherPoint()
  {
    std::vector<Point> vals;
    Real myrank = TestCommWorld->rank();
    TestCommWorld->allgather(Point(myrank, myrank+0.25, myrank+0.5),vals);

    const std::size_t comm_size = TestCommWorld->size();
    const std::size_t vec_size  = vals.size();
    CPPUNIT_ASSERT_EQUAL( comm_size, vec_size );
    for (processor_id_type i=0; i<vals.size(); i++)
      {
        Real theirrank = i;
        CPPUNIT_ASSERT_EQUAL( theirrank,            vals[i](0) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.25), vals[i](1) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.5),  vals[i](2) );
      }
  }



  void testAllGatherPairPointPoint()
  {
    std::vector<std::pair<Point, Point> > vals;
    Real myrank = TestCommWorld->rank();
    TestCommWorld->allgather
      (std::make_pair(Point(myrank, myrank+0.125, myrank+0.25), Point(myrank+0.5, myrank+0.625, myrank+0.75)), vals);

    const std::size_t comm_size = TestCommWorld->size();
    const std::size_t vec_size  = vals.size();
    CPPUNIT_ASSERT_EQUAL( comm_size, vec_size );

    for (processor_id_type i=0; i<vals.size(); i++)
      {
        Real theirrank = i;
        CPPUNIT_ASSERT_EQUAL( theirrank,             vals[i].first(0) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.125), vals[i].first(1) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.25),  vals[i].first(2) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.5),   vals[i].second(0) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.625), vals[i].second(1) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.75),  vals[i].second(2) );
      }
  }


  void testAllGatherPairRealPoint()
  {
    std::vector<std::pair<Real, Point> > vals;
    Real myrank = TestCommWorld->rank();
    TestCommWorld->allgather
      (std::make_pair(Real(myrank+0.75), Point(myrank, myrank+0.25, myrank+0.5)), vals);

    const std::size_t comm_size = TestCommWorld->size();
    const std::size_t vec_size  = vals.size();
    CPPUNIT_ASSERT_EQUAL( comm_size, vec_size );

    for (processor_id_type i=0; i<vals.size(); i++)
      {
        Real theirrank = i;
        CPPUNIT_ASSERT_EQUAL( theirrank,            vals[i].second(0) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.25), vals[i].second(1) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.5),  vals[i].second(2) );
        CPPUNIT_ASSERT_EQUAL( theirrank+Real(0.75), vals[i].first );
      }
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

        for (std::size_t i=0; i<src_val.size(); i++)
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

        for (std::size_t i=0; i<src_val.size(); i++)
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

        for (std::size_t i=0; i<src_val.size(); i++)
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

        for (std::size_t i=0; i<src_val.size(); i++)
          CPPUNIT_ASSERT_EQUAL( src_val[i] , recv_val[i] );

        // Restore default communication
        TestCommWorld->send_mode(Parallel::Communicator::DEFAULT);
      }
  }




};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelPointTest );
