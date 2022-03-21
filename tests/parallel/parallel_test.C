#include <libmesh/parallel.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class ParallelTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ParallelTest );

  CPPUNIT_TEST( testGather );
  CPPUNIT_TEST( testAllGather );
  CPPUNIT_TEST( testGatherString );
  CPPUNIT_TEST( testAllGatherString );
  CPPUNIT_TEST( testAllGatherVectorString );
  CPPUNIT_TEST( testAllGatherEmptyVectorString );
  CPPUNIT_TEST( testAllGatherHalfEmptyVectorString );
  CPPUNIT_TEST( testBroadcast );
  CPPUNIT_TEST( testBroadcastNestedType );
  CPPUNIT_TEST( testScatter );
  CPPUNIT_TEST( testBarrier );
  CPPUNIT_TEST( testMin );
  CPPUNIT_TEST( testMax );
  CPPUNIT_TEST( testMinloc );
  CPPUNIT_TEST( testMaxloc );
  CPPUNIT_TEST( testMinlocReal );
  CPPUNIT_TEST( testMaxlocReal );
  CPPUNIT_TEST( testInfinityMin );
  CPPUNIT_TEST( testInfinityMax );
  CPPUNIT_TEST( testIsendRecv );
  CPPUNIT_TEST( testIrecvSend );
  CPPUNIT_TEST( testRecvIsendSets );
  CPPUNIT_TEST( testRecvIsendVecVecs );
  CPPUNIT_TEST( testSendRecvVecVecs );
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
    LOG_UNIT_TEST;

    std::vector<processor_id_type> vals;
    TestCommWorld->gather(0,cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



  void testGatherString()
  {
    LOG_UNIT_TEST;

    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + _number[TestCommWorld->rank() % 10], vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] , vals[i] );
  }



  void testAllGather()
  {
    LOG_UNIT_TEST;

    std::vector<processor_id_type> vals;
    TestCommWorld->allgather(cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      CPPUNIT_ASSERT_EQUAL( i , vals[i] );
  }



  void testAllGatherString()
  {
    LOG_UNIT_TEST;

    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + _number[TestCommWorld->rank() % 10], vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] , vals[i] );
  }



  void testAllGatherVectorString()
  {
    LOG_UNIT_TEST;

    std::vector<std::string> vals;
    vals.push_back("Processor" + _number[TestCommWorld->rank() % 10] + "A");
    vals.push_back("Processor" + _number[TestCommWorld->rank() % 10] + "B");
    TestCommWorld->allgather(vals);

    for (processor_id_type i=0; i<(vals.size()/2); i++)
      {
        CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] + "A" , vals[2*i] );
        CPPUNIT_ASSERT_EQUAL( "Processor" + _number[i % 10] + "B" , vals[2*i+1] );
      }
  }



  void testAllGatherEmptyVectorString()
  {
    LOG_UNIT_TEST;

    std::vector<std::string> vals;
    TestCommWorld->allgather(vals);

    CPPUNIT_ASSERT( vals.empty() );
  }



  void testAllGatherHalfEmptyVectorString()
  {
    LOG_UNIT_TEST;

    std::vector<std::string> vals;

    if (!TestCommWorld->rank())
      vals.push_back("Proc 0 only");

    TestCommWorld->allgather(vals);

    CPPUNIT_ASSERT_EQUAL( vals[0], std::string("Proc 0 only") );
  }



  void testBroadcast()
  {
    LOG_UNIT_TEST;

    // Workaround for spurious warning from operator=
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100366
    std::vector<unsigned int> src{0,1,2}, dest(3,0);

    if (TestCommWorld->rank() == 0)
      dest = src;

    TestCommWorld->broadcast(dest);

    for (std::size_t i=0; i<3; i++)
      CPPUNIT_ASSERT_EQUAL( src[i] , dest[i] );
  }



  void testBroadcastNestedType()
  {
    LOG_UNIT_TEST;

    using std::pair;
    typedef pair<pair<pair<pair<int, int>, int>, int>, int> pppp;
    std::vector<pppp> src(3), dest(3);

    src[0].first.first.first.first=0;
    src[0].first.first.first.second=-1;
    src[0].first.second = -2;
    src[0].second = -3;
    src[1].first.first.first.first=10;
    src[1].first.first.first.second=9;
    src[1].first.second = 8;
    src[1].second = 7;
    src[2].first.first.first.first=20;
    src[2].first.first.first.second=19;
    src[2].first.second = 18;
    src[2].second = 17;

    if (TestCommWorld->rank() == 0)
      dest = src;

    TestCommWorld->broadcast(dest);

    for (std::size_t i=0; i<src.size(); i++)
      {
        CPPUNIT_ASSERT_EQUAL(src[i].first.first.first.first,
                             dest[i].first.first.first.first);
        CPPUNIT_ASSERT_EQUAL(src[i].first.first.first.second,
                             dest[i].first.first.first.second);
        CPPUNIT_ASSERT_EQUAL(src[i].first.first.second,
                             dest[i].first.first.second);
        CPPUNIT_ASSERT_EQUAL(src[i].first.second,
                             dest[i].first.second);
        CPPUNIT_ASSERT_EQUAL(src[i].second,
                             dest[i].second);
      }
  }



  void testScatter()
  {
    LOG_UNIT_TEST;

    // Test Scalar scatter
    {
      std::vector<processor_id_type> src;
      processor_id_type dest;

      if (TestCommWorld->rank() == 0)
        {
          src.resize(TestCommWorld->size());
          for (processor_id_type i=0; i<src.size(); i++)
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
      std::vector<std::vector<unsigned int>> src;
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
    LOG_UNIT_TEST;

    TestCommWorld->barrier();
  }



  void testMin ()
  {
    LOG_UNIT_TEST;

    unsigned int min = TestCommWorld->rank();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, static_cast<unsigned int>(0));
  }



  void testMax ()
  {
    LOG_UNIT_TEST;

    processor_id_type max = TestCommWorld->rank();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (cast_int<processor_id_type>(max+1),
                          cast_int<processor_id_type>(TestCommWorld->size()));
  }



  void testMinloc ()
  {
    LOG_UNIT_TEST;

    int min = (TestCommWorld->rank() + 1) % TestCommWorld->size();
    unsigned int minid = 0;

    TestCommWorld->minloc(min, minid);

    CPPUNIT_ASSERT_EQUAL (min, static_cast<int>(0));
    CPPUNIT_ASSERT_EQUAL (minid, static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMaxloc ()
  {
    LOG_UNIT_TEST;

    int max = TestCommWorld->rank();
    unsigned int maxid = 0;

    TestCommWorld->maxloc(max, maxid);

    CPPUNIT_ASSERT_EQUAL (max+1,
                          cast_int<int>(TestCommWorld->size()));
    CPPUNIT_ASSERT_EQUAL (maxid, static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMinlocReal ()
  {
    LOG_UNIT_TEST;

    Real min = (TestCommWorld->rank() + 1) % TestCommWorld->size();
    unsigned int minid = 0;

    TestCommWorld->minloc(min, minid);

    CPPUNIT_ASSERT_EQUAL (min, Real(0));
    CPPUNIT_ASSERT_EQUAL (minid, static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMaxlocReal ()
  {
    LOG_UNIT_TEST;

    Real max = TestCommWorld->rank();
    unsigned int maxid = 0;

    TestCommWorld->maxloc(max, maxid);

    // Hope nobody uses 1677216 procs with single precision
    CPPUNIT_ASSERT_EQUAL (max+1, Real(TestCommWorld->size()));
    CPPUNIT_ASSERT_EQUAL (maxid, static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testInfinityMin ()
  {
    LOG_UNIT_TEST;

    double min = std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, std::numeric_limits<double>::infinity());

    min = -std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    CPPUNIT_ASSERT_EQUAL (min, -std::numeric_limits<double>::infinity());
  }



  void testInfinityMax ()
  {
    LOG_UNIT_TEST;

    double max = std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (max, std::numeric_limits<double>::infinity());

    max = -std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    CPPUNIT_ASSERT_EQUAL (max, -std::numeric_limits<double>::infinity());
  }



  void testIsendRecv ()
  {
    LOG_UNIT_TEST;

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
    LOG_UNIT_TEST;

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
    LOG_UNIT_TEST;

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



  void testRecvIsendVecVecs ()
  {
    LOG_UNIT_TEST;

    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    std::vector<std::vector<unsigned int> > src_val(3), recv_val;

    src_val[0].push_back(4);  // Chosen by fair dice roll
    src_val[2].push_back(procup);
    src_val[2].push_back(TestCommWorld->rank());

    Parallel::Request request;

    if (TestCommWorld->size() > 1)
      {
        TestCommWorld->send (procup, src_val, request);

        TestCommWorld->receive (procdown,
                                recv_val);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (std::size_t i = 0; i != 3; ++i)
          CPPUNIT_ASSERT_EQUAL ( src_val[i].size(), recv_val[i].size() );

        CPPUNIT_ASSERT_EQUAL ( recv_val[0][0], static_cast<unsigned int> (4) );
        CPPUNIT_ASSERT_EQUAL ( recv_val[2][0], static_cast<unsigned int> (TestCommWorld->rank()) );
        CPPUNIT_ASSERT_EQUAL ( recv_val[2][1], procdown );

        Parallel::wait (request);

        recv_val.clear();
      }
  }


  void testSendRecvVecVecs ()
  {
    LOG_UNIT_TEST;

    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    // Any odd processor out does nothing
    if ((TestCommWorld->size() % 2) && procup == 0)
      return;

    std::vector<std::vector<unsigned int> > src_val(3), recv_val;

    src_val[0].push_back(4);  // Chosen by fair dice roll
    src_val[2].push_back(procup);
    src_val[2].push_back(TestCommWorld->rank());

    // Other even numbered processors send
    if (TestCommWorld->rank() % 2 == 0)
      TestCommWorld->send (procup, src_val);
    // Other odd numbered processors receive
    else
      {
        TestCommWorld->receive (procdown,
                                recv_val);

        CPPUNIT_ASSERT_EQUAL ( src_val.size() , recv_val.size() );

        for (std::size_t i = 0; i != 3; ++i)
          CPPUNIT_ASSERT_EQUAL ( src_val[i].size(), recv_val[i].size() );

        CPPUNIT_ASSERT_EQUAL ( recv_val[0][0], static_cast<unsigned int> (4) );
        CPPUNIT_ASSERT_EQUAL ( recv_val[2][0], static_cast<unsigned int> (TestCommWorld->rank()) );
        CPPUNIT_ASSERT_EQUAL ( recv_val[2][1], procdown );

        recv_val.clear();
      }
  }



  void testSemiVerify ()
  {
    LOG_UNIT_TEST;

    double inf = std::numeric_limits<double>::infinity();

    double *infptr = TestCommWorld->rank()%2 ? NULL : &inf;

    CPPUNIT_ASSERT (TestCommWorld->semiverify(infptr));

    inf = -std::numeric_limits<double>::infinity();

    CPPUNIT_ASSERT (TestCommWorld->semiverify(infptr));
  }


  void testSplit ()
  {
    LOG_UNIT_TEST;

    Parallel::Communicator subcomm;
    unsigned int rank = TestCommWorld->rank();
    unsigned int color = rank % 2;
    TestCommWorld->split(color, rank, subcomm);

    CPPUNIT_ASSERT(subcomm.size() >= 1);
    CPPUNIT_ASSERT(subcomm.size() >= TestCommWorld->size() / 2);
    CPPUNIT_ASSERT(subcomm.size() <= TestCommWorld->size() / 2 + 1);
  }


  void testSplitByType ()
  {
    LOG_UNIT_TEST;

    Parallel::Communicator subcomm;
    unsigned int rank = TestCommWorld->rank();
    Parallel::info i = 0;
    int type = 0;
#ifdef LIBMESH_HAVE_MPI
    type = MPI_COMM_TYPE_SHARED;
    i = MPI_INFO_NULL;
#endif
    TestCommWorld->split_by_type(type, rank, i, subcomm);

    CPPUNIT_ASSERT(subcomm.size() >= 1);
    CPPUNIT_ASSERT(subcomm.size() <= TestCommWorld->size());
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelTest );
