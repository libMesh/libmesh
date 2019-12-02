#include <timpi/timpi.h>

#define TIMPI_UNIT_ASSERT(expr) \
  if (!(expr)) \
    timpi_error();

using namespace TIMPI;

Communicator *TestCommWorld;

std::vector<std::string> pt_number;

  void setUp()
  {
    pt_number.resize(10);
    pt_number[0] = "Zero";
    pt_number[1] = "One";
    pt_number[2] = "Two";
    pt_number[3] = "Three";
    pt_number[4] = "Four";
    pt_number[5] = "Five";
    pt_number[6] = "Six";
    pt_number[7] = "Seven";
    pt_number[8] = "Eight";
    pt_number[9] = "Nine";
  }


  void testGather()
  {
    std::vector<processor_id_type> vals;
    TestCommWorld->gather(0,cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        TIMPI_UNIT_ASSERT( i  == vals[i] );
  }



  void testGatherString()
  {
    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + pt_number[TestCommWorld->rank() % 10], vals);

    if (TestCommWorld->rank() == 0)
      for (processor_id_type i=0; i<vals.size(); i++)
        TIMPI_UNIT_ASSERT( "Processor" + pt_number[i % 10]  == vals[i] );
  }



  void testAllGather()
  {
    std::vector<processor_id_type> vals;
    TestCommWorld->allgather(cast_int<processor_id_type>(TestCommWorld->rank()),vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      TIMPI_UNIT_ASSERT( i  == vals[i] );
  }



  void testAllGatherString()
  {
    std::vector<std::string> vals;
    TestCommWorld->gather(0, "Processor" + pt_number[TestCommWorld->rank() % 10], vals);

    for (processor_id_type i=0; i<vals.size(); i++)
      TIMPI_UNIT_ASSERT( "Processor" + pt_number[i % 10]  == vals[i] );
  }



  void testAllGatherVectorString()
  {
    std::vector<std::string> vals;
    vals.push_back("Processor" + pt_number[TestCommWorld->rank() % 10] + "A");
    vals.push_back("Processor" + pt_number[TestCommWorld->rank() % 10] + "B");
    TestCommWorld->allgather(vals);

    for (processor_id_type i=0; i<(vals.size()/2); i++)
      {
        TIMPI_UNIT_ASSERT( "Processor" + pt_number[i % 10] + "A"  == vals[2*i] );
        TIMPI_UNIT_ASSERT( "Processor" + pt_number[i % 10] + "B"  == vals[2*i+1] );
      }
  }



  void testAllGatherEmptyVectorString()
  {
    std::vector<std::string> vals;
    TestCommWorld->allgather(vals);

    TIMPI_UNIT_ASSERT( vals.empty() );
  }



  void testAllGatherHalfEmptyVectorString()
  {
    std::vector<std::string> vals;

    if (!TestCommWorld->rank())
      vals.push_back("Proc 0 only");

    TestCommWorld->allgather(vals);

    TIMPI_UNIT_ASSERT( vals[0] == std::string("Proc 0 only") );
  }



  void testBroadcast(const bool use_proc_null)
  {
    std::vector<unsigned int> src(3), dest(3);

    src[0]=0;
    src[1]=1;
    src[2]=2;

    if (TestCommWorld->rank() == 0)
      dest = src;

    unsigned int root_id = 0;
#ifdef LIBMESH_HAVE_MPI
    if (use_proc_null)
      root_id = (unsigned int)(TestCommWorld->rank() == 0 ? MPI_ROOT : MPI_PROC_NULL);
#endif

    TestCommWorld->broadcast(dest, root_id);

    for (std::size_t i=0; i<src.size(); i++)
      TIMPI_UNIT_ASSERT( src[i]  == dest[i] );
  }



  void testBroadcastNestedType(const bool use_proc_null)
  {
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

    unsigned int root_id = 0;
  #ifdef LIBMESH_HAVE_MPI
      if (use_proc_null)
        root_id = (unsigned int)(TestCommWorld->rank() == 0 ? MPI_ROOT : MPI_PROC_NULL);
  #endif

    TestCommWorld->broadcast(dest, root_id);

    for (std::size_t i=0; i<src.size(); i++)
      {
        TIMPI_UNIT_ASSERT(src[i].first.first.first.first ==
                             dest[i].first.first.first.first);
        TIMPI_UNIT_ASSERT(src[i].first.first.first.second ==
                             dest[i].first.first.first.second);
        TIMPI_UNIT_ASSERT(src[i].first.first.second ==
                             dest[i].first.first.second);
        TIMPI_UNIT_ASSERT(src[i].first.second ==
                             dest[i].first.second);
        TIMPI_UNIT_ASSERT(src[i].second ==
                             dest[i].second);
      }
  }



  void testScatter()
  {
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

      TIMPI_UNIT_ASSERT( TestCommWorld->rank() == dest );
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
        TIMPI_UNIT_ASSERT( TestCommWorld->rank() * CHUNK_SIZE + i == dest[i] );
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
        TIMPI_UNIT_ASSERT( start_value + i == dest[i] );
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
        TIMPI_UNIT_ASSERT( start_value + i == dest[i] );
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

    TIMPI_UNIT_ASSERT (min == static_cast<unsigned int>(0));
  }



  void testMax ()
  {
    processor_id_type max = TestCommWorld->rank();

    TestCommWorld->max(max);

    TIMPI_UNIT_ASSERT (cast_int<processor_id_type>(max+1) ==
                          cast_int<processor_id_type>(TestCommWorld->size()));
  }



  void testMinloc ()
  {
    int min = (TestCommWorld->rank() + 1) % TestCommWorld->size();
    unsigned int minid = 0;

    TestCommWorld->minloc(min, minid);

    TIMPI_UNIT_ASSERT (min == static_cast<int>(0));
    TIMPI_UNIT_ASSERT (minid == static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMaxloc ()
  {
    int max = TestCommWorld->rank();
    unsigned int maxid = 0;

    TestCommWorld->maxloc(max, maxid);

    TIMPI_UNIT_ASSERT (max+1 ==
                          cast_int<int>(TestCommWorld->size()));
    TIMPI_UNIT_ASSERT (maxid == static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMinlocDouble ()
  {
    double min = (TestCommWorld->rank() + 1) % TestCommWorld->size();
    unsigned int minid = 0;

    TestCommWorld->minloc(min, minid);

    TIMPI_UNIT_ASSERT (min == double(0));
    TIMPI_UNIT_ASSERT (minid == static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testMaxlocDouble ()
  {
    double max = TestCommWorld->rank();
    unsigned int maxid = 0;

    TestCommWorld->maxloc(max, maxid);

    // Hope nobody uses 1677216 procs with single precision
    TIMPI_UNIT_ASSERT (max+1 == double(TestCommWorld->size()));
    TIMPI_UNIT_ASSERT (maxid == static_cast<unsigned int>(TestCommWorld->size()-1));
  }



  void testInfinityMin ()
  {
    double min = std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    TIMPI_UNIT_ASSERT (min == std::numeric_limits<double>::infinity());

    min = -std::numeric_limits<double>::infinity();

    TestCommWorld->min(min);

    TIMPI_UNIT_ASSERT (min == -std::numeric_limits<double>::infinity());
  }



  void testInfinityMax ()
  {
    double max = std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    TIMPI_UNIT_ASSERT (max == std::numeric_limits<double>::infinity());

    max = -std::numeric_limits<double>::infinity();

    TestCommWorld->max(max);

    TIMPI_UNIT_ASSERT (max == -std::numeric_limits<double>::infinity());
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

    TIMPI::Request request;

    if (TestCommWorld->size() > 1)
      {
        // Default communication
        TestCommWorld->send_mode(TIMPI::Communicator::DEFAULT);

        TestCommWorld->send (procup,
                             src_val,
                             request);

        TestCommWorld->receive (procdown,
                                recv_val);

        TIMPI::wait (request);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i=0; i<src_val.size(); i++)
          TIMPI_UNIT_ASSERT( src_val[i]  == recv_val[i] );


        // Synchronous communication
        TestCommWorld->send_mode(TIMPI::Communicator::SYNCHRONOUS);
        std::fill (recv_val.begin(), recv_val.end(), 0);

        TestCommWorld->send (procup,
                             src_val,
                             request);

        TestCommWorld->receive (procdown,
                                recv_val);

        TIMPI::wait (request);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i=0; i<src_val.size(); i++)
          TIMPI_UNIT_ASSERT( src_val[i]  == recv_val[i] );

        // Restore default communication
        TestCommWorld->send_mode(TIMPI::Communicator::DEFAULT);
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

    TIMPI::Request request;

    if (TestCommWorld->size() > 1)
      {
        // Default communication
        TestCommWorld->send_mode(TIMPI::Communicator::DEFAULT);

        TestCommWorld->receive (procdown,
                                recv_val,
                                request);

        TestCommWorld->send (procup,
                             src_val);

        TIMPI::wait (request);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i=0; i<src_val.size(); i++)
          TIMPI_UNIT_ASSERT( src_val[i]  == recv_val[i] );

        // Synchronous communication
        TestCommWorld->send_mode(TIMPI::Communicator::SYNCHRONOUS);
        std::fill (recv_val.begin(), recv_val.end(), 0);


        TestCommWorld->receive (procdown,
                                recv_val,
                                request);

        TestCommWorld->send (procup,
                             src_val);

        TIMPI::wait (request);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i=0; i<src_val.size(); i++)
          TIMPI_UNIT_ASSERT( src_val[i]  == recv_val[i] );

        // Restore default communication
        TestCommWorld->send_mode(TIMPI::Communicator::DEFAULT);
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

    TIMPI::Request request;

    if (TestCommWorld->size() > 1)
      {
        TestCommWorld->send (procup, src_val, request);

        TestCommWorld->receive (procdown,
                                recv_val);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::set<unsigned int>::const_iterator
               it = src_val.begin(), end = src_val.end(); it != end;
             ++it)
          TIMPI_UNIT_ASSERT ( recv_val.count(*it) );

        TIMPI::wait (request);

        recv_val.clear();
      }
  }



  void testRecvIsendVecVecs ()
  {
    unsigned int procup = (TestCommWorld->rank() + 1) %
      TestCommWorld->size();
    unsigned int procdown = (TestCommWorld->size() +
                             TestCommWorld->rank() - 1) %
      TestCommWorld->size();

    std::vector<std::vector<unsigned int> > src_val(3), recv_val;

    src_val[0].push_back(4);  // Chosen by fair dice roll
    src_val[2].push_back(procup);
    src_val[2].push_back(TestCommWorld->rank());

    TIMPI::Request request;

    if (TestCommWorld->size() > 1)
      {
        TestCommWorld->send (procup, src_val, request);

        TestCommWorld->receive (procdown,
                                recv_val);

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i = 0; i != 3; ++i)
          TIMPI_UNIT_ASSERT ( src_val[i].size() == recv_val[i].size() );

        TIMPI_UNIT_ASSERT ( recv_val[0][0] == static_cast<unsigned int> (4) );
        TIMPI_UNIT_ASSERT ( recv_val[2][0] == static_cast<unsigned int> (TestCommWorld->rank()) );
        TIMPI_UNIT_ASSERT ( recv_val[2][1] == procdown );

        TIMPI::wait (request);

        recv_val.clear();
      }
  }


  void testSendRecvVecVecs ()
  {
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

        TIMPI_UNIT_ASSERT ( src_val.size()  == recv_val.size() );

        for (std::size_t i = 0; i != 3; ++i)
          TIMPI_UNIT_ASSERT ( src_val[i].size() == recv_val[i].size() );

        TIMPI_UNIT_ASSERT ( recv_val[0][0] == static_cast<unsigned int> (4) );
        TIMPI_UNIT_ASSERT ( recv_val[2][0] == static_cast<unsigned int> (TestCommWorld->rank()) );
        TIMPI_UNIT_ASSERT ( recv_val[2][1] == procdown );

        recv_val.clear();
      }
  }



  void testSemiVerify ()
  {
    double inf = std::numeric_limits<double>::infinity();

    double *infptr = TestCommWorld->rank()%2 ? NULL : &inf;

    TIMPI_UNIT_ASSERT (TestCommWorld->semiverify(infptr));

    inf = -std::numeric_limits<double>::infinity();

    TIMPI_UNIT_ASSERT (TestCommWorld->semiverify(infptr));
  }


  void testSplit ()
  {
    TIMPI::Communicator subcomm;
    unsigned int rank = TestCommWorld->rank();
    unsigned int color = rank % 2;
    TestCommWorld->split(color, rank, subcomm);

    TIMPI_UNIT_ASSERT(subcomm.size() >= 1);
    TIMPI_UNIT_ASSERT(subcomm.size() >= TestCommWorld->size() / 2);
    TIMPI_UNIT_ASSERT(subcomm.size() <= TestCommWorld->size() / 2 + 1);
  }


  void testSplitByType ()
  {
    TIMPI::Communicator subcomm;
    unsigned int rank = TestCommWorld->rank();
    TIMPI::info i = 0;
    int type = 0;
#ifdef LIBMESH_HAVE_MPI
    type = MPI_COMM_TYPE_SHARED;
    i = MPI_INFO_NULL;
#endif
    TestCommWorld->split_by_type(type, rank, i, subcomm);

    TIMPI_UNIT_ASSERT(subcomm.size() >= 1);
    TIMPI_UNIT_ASSERT(subcomm.size() <= TestCommWorld->size());
  }


int main(int argc, const char * const * argv)
{
  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  setUp();

  testGather();
  testAllGather();
  testGatherString();
  testAllGatherString();
  testAllGatherVectorString();
  testAllGatherEmptyVectorString();
  testAllGatherHalfEmptyVectorString();
  testBroadcast(/* use_proc_null = */ false);
  testBroadcast(/* use_proc_null = */ true);
  testBroadcastNestedType(/* use_proc_null = */ false);
  testBroadcastNestedType(/* use_proc_null = */ true);
  testScatter();
  testBarrier();
  testMin();
  testMax();
  testMinloc();
  testMaxloc();
  testMinlocDouble();
  testMaxlocDouble();
  testInfinityMin();
  testInfinityMax();
  testIsendRecv();
  testIrecvSend();
  testRecvIsendSets();
  testRecvIsendVecVecs();
  testSendRecvVecVecs();
  testSemiVerify();
  testSplit();

  return 0;
}
