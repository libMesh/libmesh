// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include "libmesh/mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/checkpoint_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/metis_partitioner.h"

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

class CheckpointIOTest : public CppUnit::TestCase {
  /**
   * This test verifies that we can write files with the CheckpointIO object.
   */
public:
  CPPUNIT_TEST_SUITE( CheckpointIOTest );

  CPPUNIT_TEST( testSplitter );

  CPPUNIT_TEST_SUITE_END();

protected:

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  // Test that we can write multiple checkpoint files from a single processor.
  void testSplitter()
  {
    // In this test, we partition the mesh into n_procs parts.
    const unsigned int n_procs = 2;

    // The number of elements in the original mesh.  For verification
    // later.
    dof_id_type original_n_elem = 0;

    {
      MetisPartitioner partitioner;
      ReplicatedMesh mesh(*TestCommWorld);

      MeshTools::Generation::build_square(mesh,
                                          4,  4,
                                          0., 1.,
                                          0., 1.,
                                          QUAD4);

      // Store the number of elements that were in the original mesh.
      original_n_elem = mesh.n_elem();

      // Partition the mesh into n_procs pieces
      partitioner.partition(mesh, n_procs);

      // Write out checkpoint files for each processor.
      for (unsigned i=0; i<n_procs; ++i)
        {
          CheckpointIO cpr(mesh);
          cpr.current_processor_id() = i;
          cpr.current_n_processors() = n_procs;
          cpr.binary() = true;
          cpr.parallel() = true;
          cpr.write("checkpoint_splitter.cpr");
        }
    }

    // Test that we can read in the files we wrote and sum up to the
    // same total number of elements.
    {
      unsigned int read_in_elements = 0;

      for (unsigned i=0; i<n_procs; ++i)
        {
          Mesh mesh(*TestCommWorld);
          CheckpointIO cpr(mesh);
          cpr.current_processor_id() = i;
          cpr.current_n_processors() = n_procs;
          cpr.binary() = true;
          cpr.parallel() = true;
          cpr.read("checkpoint_splitter.cpr");
          read_in_elements += std::distance(mesh.pid_elements_begin(i),
                                            mesh.pid_elements_end(i));
        }

      // Verify that we read in exactly as many elements on each proc as we started with.
      CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(read_in_elements), original_n_elem);
    }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( CheckpointIOTest );
