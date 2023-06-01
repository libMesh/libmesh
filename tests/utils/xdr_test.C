// Unit test includes
#include "libmesh_cppunit.h"
#include "test_comm.h"

// libMesh includes
#include <libmesh/libmesh.h>
#include <libmesh/xdr_cxx.h>
#include <libmesh/int_range.h>
#include <timpi/communicator.h>

// C++ includes
#include <vector>
#include <sstream>

using namespace libMesh;

class XdrTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( XdrTest );

  CPPUNIT_TEST( testDataVec );
  CPPUNIT_TEST( testDataStream );

  CPPUNIT_TEST_SUITE_END();

private:
  template <typename ReadLambda, typename WriteLambda>
  void test_read_write(ReadLambda & act_read, WriteLambda & act_write)
  {
    // There was once a weird bug mutating our TestCommWorld in
    // a previous unit test, which turned *this* test into a race
    // condition.  Let's make sure that never happens again.
    CPPUNIT_ASSERT_EQUAL(TestCommWorld->rank(),
                         libMesh::global_processor_id());

    // Test reading/writing on 1 processor
    if (TestCommWorld->rank() == 0)
      {
        // Write to file
        {
          Xdr xdr("output.dat", WRITE);
          act_write(xdr);
        }

        // Read from file
        {
          Xdr xdr("output.dat", READ);
          act_read(xdr);
        }
      }

    // Test reading/writing to a stream directly on all processors
    {
      std::stringstream stream;

      // Write to stream
      {
        std::ostream ostream(stream.rdbuf());
        Xdr xdr(ostream);
        act_write(xdr);
      }

      // Rewind stream
      stream.seekp(0);

      // Read from stream
      {
        std::istream istream(stream.rdbuf());
        Xdr xdr(istream);
        act_read(xdr);
      }
    }
  }

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testDataVec ()
  {
    LOG_UNIT_TEST;

    std::vector<Real> vec = {libMesh::pi, 2*libMesh::pi, 3*libMesh::pi};

    auto act_write = [&vec](Xdr & xdr) { xdr.data(vec, "# This is a comment"); };
    auto act_read = [&vec](Xdr & xdr) {
      std::vector<Real> vec_in;
      xdr.data(vec_in);

      // Check that correct number of values were read in
      CPPUNIT_ASSERT_EQUAL(vec_in.size(), vec.size());

      // Check that values were written/read with sufficient accuracy
      for (auto i : index_range(vec_in))
        LIBMESH_ASSERT_FP_EQUAL(vec[i], vec_in[i], TOLERANCE); };

    test_read_write(act_read, act_write);
  }

  void testDataStream ()
  {
    LOG_UNIT_TEST;

    std::vector<Real> vec(100);
    for (auto i : index_range(vec))
      vec[i] = static_cast<Real>(i+1) / vec.size();

    // Write. If "line_break" does not exactly divide the vector size,
    // there will be one line with fewer values than the others.
    auto act_write = [&vec](Xdr & xdr) {
      xdr.data_stream(vec.data(), vec.size(), /*line_break=*/16); };
    // Read. To use data_stream(), the storage needs to be pre-sized
    // and the line_break parameter is ignored.
    auto act_read = [&vec](Xdr & xdr) {
      std::vector<Real> vec_in(100);
      xdr.data_stream(vec_in.data(), vec_in.size());

      // Check that values were written/read correctly
      for (auto i : index_range(vec_in))
        LIBMESH_ASSERT_FP_EQUAL(vec[i], vec_in[i], TOLERANCE); };

    test_read_write(act_read, act_write);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( XdrTest );
