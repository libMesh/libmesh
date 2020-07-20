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

using namespace libMesh;

class XdrTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( XdrTest );

  CPPUNIT_TEST( testDataVec );
  CPPUNIT_TEST( testDataStream );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testDataVec ()
  {
    std::vector<Real> vec = {libMesh::pi, 2*libMesh::pi, 3*libMesh::pi};

    // Test reading/writing on 1 processor
    if (TestCommWorld->rank() == 0)
      {
        // Write to file
        {
          Xdr xdr("output.dat", WRITE);
          xdr.data(vec, "# This is a comment");
        }

        // Read from file
        {
          Xdr xdr("output.dat", READ);
          std::vector<Real> vec_in;
          xdr.data(vec_in);

          // Check that correct number of values were read in
          CPPUNIT_ASSERT_EQUAL(vec_in.size(), vec.size());

          // Check that values were written/read with sufficient accuracy
          for (auto i : index_range(vec_in))
            LIBMESH_ASSERT_FP_EQUAL(vec[i], vec_in[i], TOLERANCE);
        }
      }
  }

  void testDataStream ()
  {
    std::vector<Real> vec(100);
    for (auto i : index_range(vec))
      vec[i] = static_cast<Real>(i+1) / vec.size();

    // Test reading/writing on 1 processor
    if (TestCommWorld->rank() == 0)
      {
        // Write to file. If "line_break" does not exactly divide the
        // vector size, there will be one line with fewer values than
        // the others.
        {
          Xdr xdr("output.dat", WRITE);
          xdr.data_stream(vec.data(), vec.size(), /*line_break=*/16);
        }

        // Read from file. To use data_stream(), the storage needs to
        // be pre-sized and the line_break parameter is ignored.
        {
          Xdr xdr("output.dat", READ);
          std::vector<Real> vec_in(100);
          xdr.data_stream(vec_in.data(), vec_in.size());

          // Check that values were written/read correctly
          for (auto i : index_range(vec_in))
            LIBMESH_ASSERT_FP_EQUAL(vec[i], vec_in[i], TOLERANCE);
        }
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( XdrTest );
