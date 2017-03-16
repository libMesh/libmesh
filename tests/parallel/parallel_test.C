// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel.h>

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

class ParallelTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParallelTest );

  CPPUNIT_TEST( testGather );
  CPPUNIT_TEST( testAllGather );
  CPPUNIT_TEST( testGatherString );
  CPPUNIT_TEST( testAllGatherString );
  CPPUNIT_TEST( testBroadcast );
  CPPUNIT_TEST( testScatter );
  CPPUNIT_TEST( testBarrier );
  CPPUNIT_TEST( testMin );
  CPPUNIT_TEST( testMax );
  CPPUNIT_TEST( testInfinityMin );
  CPPUNIT_TEST( testInfinityMax );
  CPPUNIT_TEST( testIsendRecv );
  CPPUNIT_TEST( testIrecvSend );
  CPPUNIT_TEST( testRecvIsendSets );
  CPPUNIT_TEST( testSemiVerify );
  CPPUNIT_TEST( testSplit );

  CPPUNIT_TEST_SUITE_END();

private:
  std::vector<std::string> _number;

public:
  void setUp()
  {
    _number.resize(10);
    _number[0] = "Zero";
    _number[1] = "One";
    _number[2] = "Two";
    _number[3] = "Three";
    _number[4] = "Four";
    _number[5] = "Five";
    _number[6] = "Six";
    _number[7] = "Seven";
    _number[8] = "Eight";
    _number[9] = "Nine";
  }

  void tearDown()
  {}



  void testGather()
  {
    std::vector<processor_id_type> vals;
    TestCommWorld->gather(0,cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



  void testGatherString()
  {
    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + _number[TestCommWorld->rank() % 10], vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] , vals[i] );
  }



  void testAllGather()
  {
    std::vector<processor_id_type> vals;
    TestCommWorld->allgather(cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



  void testAllGatherString()
  {
    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + _number[TestCommWorld->rank() % 10], vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] , vals[i] );
  }



  void testBroadcast()
  {
    std::vector<unsigned int> src(3), dest(3);

    src[0]=0;
    src[1]=1;
    src[2]=2;

    if (TestCommWorld->rank() == 0)
      dest = src;

    TestCommWorld->broadcast(dest);

    for (std::size_t i=0; i<src.size(); i++)
      CPPUNIT_ASSERT_EQUAL( src[i] , dest[i] );
  }



  void testScatter()
  {
    // Test Scalar scatter
    {
      std::vector<unsigned int> src;
      unsigned int dest;

      if (TestCommWorld->rank() == 0)
        {
          src.resize(TestCommWorld->size());
          for (std::size_t i=0; i<src.size(); i++)
            src[i] = i;
        }

      TestCommWorld->scatter(src, dest);

      CPPUNIT_ASSERT_EQUAL( TestCommWorld->rank(), dest );
    }

    // Test Vector Scatter (equal-sized chunks)
    {
      std::vector<unsigned int> src;
      std::vector<unsigned int> dest;
      static const unsigned int CHUNK_SIZE = 3;

      if (TestCommWorld->rank() == 0)
        {
          src.resize(TestCommWorld->size() * CHUNK_SIZE);
          for (std::size_t i=0; i<src.size(); i++)
            src[i] = i;
        }

      TestCommWorld->scatter(src, dest);

      for (unsigned int i=0; i<CHUNK_SIZE; i++)
        CPPUNIT_ASSERT_EQUAL( TestCommWorld->rank() * CHUNK_SIZE + i, dest[i] );
    }

    // Test Vector Scatter (jagged chunks)
    {
      std::vector<unsigned int> src;
      std::vector<unsigned int> dest;
      std::vector<int> counts;

      if (TestCommWorld->rank() == 0)
        {
          // Give each processor "rank" number of items ( Sum i=1..n == (n * (n + 1))/2 )
          src.resize((TestCommWorld->size() * (TestCommWorld->size() + 1)) / 2);
          counts.resize(TestCommWorld->size());

          for (std::size_t i=0; i<src.size(); i++)
            src[i] = i;
          for (unsigned int i=0; i<TestCommWorld->size(); i++)
            counts[i] = static_cast<int>(i+1);
        }

      TestCommWorld->scatter(src, counts, dest);

      unsigned int start_value = (TestCommWorld->rank() * (TestCommWorld->rank() + 1)) / 2;
      for (unsigned int i=0; i<=TestCommWorld->rank(); i++)
        CPPUNIT_ASSERT_EQUAL( start_value + i, dest[i] );
    }

    // Test Vector of Vector Scatter
    {
      std::vector<std::vector<unsigned int> > src;
      std::vector<unsigned int> dest;

      if (TestCommWorld->rank() == 0)
        {
          // Give each processor "rank" number of items ( Sum i=1..n == (n * (n + 1))/2 )
          src.resize(TestCommWorld->size());
          for (std::size_t i=0; i<src.size(); ++i)
            src[i].resize(i+1);

          unsigned int global_counter = 0;
          for (std::size_t i=0; i<src.size(); i++)
            for (std::size_t j=0; j<src[i].size(); j++)
              src[i][j] = global_counter++;
        }

      TestCommWorld->scatter(src, dest);

      unsigned int start_value = (TestCommWorld->rank() * (TestCommWorld->rank() + 1)) / 2;
      for (unsigned int i=0; i<=TestCommWorld->rank(); i++)
        CPPUNIT_ASSERT_EQUAL( start_value + i, dest[i] );
    }
  }



  void testBarrier()
  {
    TestCommWorld->barrier();
  }



  void testMin ()
  {
    unsigned int min = TestCommWorld->rank();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, static_cast<unsigned int>(0));
  }



  void testMax ()
  {
    processor_id_type max = TestCommWorld->rank();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (cast_int<processor_id_type>(max+1),
                          cast_int<processor_id_type>(TestCommWorld->size()));
  }



  void testInfinityMin ()
  {
    double min = std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, std::numeric_limits<double>::infinity());

    min = -std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, -std::numeric_limits<double>::infinity());
  }



  void testInfinityMax ()
  {
    double max = std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (max, std::numeric_limits<double>::infinity());

    max = -std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (max, -std::numeric_limits<double>::infinity());
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


  void testRecvIsendSets ()
  {
    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    std::set<unsigned int> src_val, recv_val;

    src_val.insert(4);  // Chosen by fair dice roll
    src_val.insert(42);
    src_val.insert(1337);

    Parallel::Request request;

    if (TestCommWorld->size() > 1)
      {
        TestCommWorld->send (procup, src_val, request);

        TestCommWorld->receive (procdown,
                                recv_val);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (std::set<unsigned int>::const_iterator
               it = src_val.begin(), end = src_val.end(); it != end;
             ++it)
          CPPUNIT_ASSERT ( recv_val.count(*it) );

        Parallel::wait (request);

        recv_val.clear();
      }
  }



  void testSemiVerify ()
  {
    double inf = std::numeric_limits<double>::infinity();

    double *infptr = TestCommWorld->rank()%2 ? NULL : &inf;

    CPPUNIT_ASSERT (TestCommWorld->semiverify(infptr));

    inf = -std::numeric_limits<double>::infinity();

    CPPUNIT_ASSERT (TestCommWorld->semiverify(infptr));
  }


  void testSplit ()
  {
    Parallel::Communicator subcomm;
    unsigned int rank = TestCommWorld->rank();
    unsigned int color = rank % 2;
    TestCommWorld->split(color, rank, subcomm);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelTest );
