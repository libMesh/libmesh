// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include "libmesh/distributed_mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/checkpoint_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/partitioner.h"

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

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testAsciiDistRepSplitter );
  CPPUNIT_TEST( testBinaryDistRepSplitter );
  CPPUNIT_TEST( testAsciiRepDistSplitter );
  CPPUNIT_TEST( testBinaryRepDistSplitter );
  CPPUNIT_TEST( testAsciiRepRepSplitter );
  CPPUNIT_TEST( testBinaryRepRepSplitter );
  CPPUNIT_TEST( testAsciiDistDistSplitter );
  CPPUNIT_TEST( testBinaryDistDistSplitter );
#endif

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
  template <typename MeshA, typename MeshB>
  void testSplitter(bool binary, bool using_distmesh)
  {
    // The CheckpointIO-based splitter requires XDR.
#ifdef LIBMESH_HAVE_XDR

    // In this test, we partition the mesh into n_procs parts.  Don't
    // try to partition a DistributedMesh into more parts than we have
    // processors, though.
    const unsigned int n_procs = using_distmesh ?
      std::min(static_cast<processor_id_type>(2), TestCommWorld->size()) :
      2;

    // The number of elements in the original mesh.  For verification
    // later.
    dof_id_type original_n_elem = 0;

    const std::string filename =
      std::string("checkpoint_splitter.cp") + (binary ? "r" : "a");

    {
      MeshA mesh(*TestCommWorld);

      MeshTools::Generation::build_square(mesh,
                                          4,  4,
                                          0., 1.,
                                          0., 1.,
                                          QUAD4);

      // Store the number of elements that were in the original mesh.
      original_n_elem = mesh.n_elem();

      // Partition the mesh into n_procs pieces
      mesh.partition(n_procs);

      // Write out checkpoint files for each piece.  Since on a
      // ReplicatedMesh we might have more pieces than we do
      // processors, some processors may have to write out more than
      // one piece.
      CheckpointIO cpr(mesh);
      cpr.current_processor_ids().clear();
      for (processor_id_type pid = mesh.processor_id(); pid < n_procs; pid += mesh.n_processors())
        cpr.current_processor_ids().push_back(pid);
      cpr.current_n_processors() = n_procs;
      cpr.binary() = binary;
      cpr.parallel() = true;
      cpr.write(filename);
    }

    TestCommWorld->barrier();

    // Test that we can read in the files we wrote and sum up to the
    // same total number of elements.
    {
      MeshB mesh(*TestCommWorld);
      CheckpointIO cpr(mesh);
      cpr.current_n_processors() = n_procs;
      cpr.binary() = binary;
      cpr.read(filename);

      std::size_t read_in_elements = 0;

      for (unsigned pid=mesh.processor_id(); pid<n_procs; pid += mesh.n_processors())
        {
          read_in_elements += std::distance(mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid));
        }
      mesh.comm().sum(read_in_elements);

      // Verify that we read in exactly as many elements as we started with.
      CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(read_in_elements), original_n_elem);
    }
#endif // LIBMESH_HAVE_XDR
  }

  void testAsciiDistRepSplitter()
  {
    testSplitter<DistributedMesh, ReplicatedMesh>(false, true);
  }

  void testBinaryDistRepSplitter()
  {
    testSplitter<DistributedMesh, ReplicatedMesh>(true, true);
  }

  void testAsciiRepDistSplitter()
  {
    testSplitter<ReplicatedMesh, DistributedMesh>(false, true);
  }

  void testBinaryRepDistSplitter()
  {
    testSplitter<ReplicatedMesh, DistributedMesh>(true, true);
  }

  void testAsciiRepRepSplitter()
  {
    testSplitter<ReplicatedMesh, ReplicatedMesh>(false, false);
  }

  void testBinaryRepRepSplitter()
  {
    testSplitter<ReplicatedMesh, ReplicatedMesh>(true, false);
  }

  void testAsciiDistDistSplitter()
  {
    testSplitter<DistributedMesh, DistributedMesh>(false, true);
  }

  void testBinaryDistDistSplitter()
  {
    testSplitter<DistributedMesh, DistributedMesh>(true, true);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( CheckpointIOTest );
