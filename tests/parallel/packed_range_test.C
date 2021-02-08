#include <libmesh/parallel.h>
#include <libmesh/null_output_iterator.h>

#include <vector>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class PackedRangeTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( PackedRangeTest );

  CPPUNIT_TEST( testNullAllGather );
  CPPUNIT_TEST( testNullSendReceive );
  CPPUNIT_TEST( testContainerSendReceive );
  //  CPPUNIT_TEST( testAdapterSendReceive );
  //  CPPUNIT_TEST( testPointerAdapterSendReceive );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



  void testNullAllGather()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1);
    if (TestCommWorld->rank() == 0)
      send[0].assign("Hello");
    else
      send[0].assign("Goodbye");

    TestCommWorld->allgather_packed_range
      ((void *)(NULL), send.begin(), send.end(),
       libMesh::null_output_iterator<std::string>());
  }


  void testNullSendReceive()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1);
    const unsigned int my_rank = TestCommWorld->rank();
    const unsigned int dest_rank =
      (my_rank + 1) % TestCommWorld->size();
    const unsigned int source_rank =
      (my_rank + TestCommWorld->size() - 1) % TestCommWorld->size();

    {
      std::ostringstream os;
      os << my_rank;
      send[0] = os.str();
    }

    TestCommWorld->send_receive_packed_range
      (dest_rank, (void *)(NULL), send.begin(), send.end(),
       source_rank, (void *)(NULL),
       libMesh::null_output_iterator<std::string>(),
       (std::string*)NULL);
  }


  void testContainerSendReceive()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1), recv;

    const unsigned int my_rank = TestCommWorld->rank();
    const unsigned int dest_rank =
      (my_rank + 1) % TestCommWorld->size();
    const unsigned int source_rank =
      (my_rank + TestCommWorld->size() - 1) % TestCommWorld->size();

    {
      std::ostringstream os;
      os << my_rank;
      send[0] = os.str();
    }

    TestCommWorld->send_receive_packed_range
      (dest_rank, (void *)(NULL), send.begin(), send.end(),
       source_rank, (void *)(NULL),
       std::back_inserter(recv),
       (std::string*)NULL);

    CPPUNIT_ASSERT_EQUAL(recv.size(), std::size_t(1));

    std::string check;
    {
      std::ostringstream os;
      os << source_rank;
      check = os.str();
    }

    CPPUNIT_ASSERT_EQUAL(recv[0], check);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PackedRangeTest );
