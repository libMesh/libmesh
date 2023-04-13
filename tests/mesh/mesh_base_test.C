
#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class MeshBaseTest : public CppUnit::TestCase {
  /**
   * This test verifies various operations on the MeshBase and derived
   * classes.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshBaseTest );

/* Tests need a 2d mesh */
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testDistributedMeshVerifyIsPrepared );
  CPPUNIT_TEST( testMeshVerifyIsPrepared );
  CPPUNIT_TEST( testReplicatedMeshVerifyIsPrepared );
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {}

  void tearDown() {}

  void testMeshBaseVerifyIsPrepared(UnstructuredMesh & mesh)
  {
    /**
     * Build a 2d 2x2 square mesh (mesh_one) covering [0.0, 1.0] x [0.0, 1.0]
     * with linear Quad elements.
     */
    MeshTools::Generation::build_square(mesh,
                                        2, 2,
                                        0., 1.,
                                        0., 1.,
                                        QUAD9);

    // build_square does its own prepare_for_use(), so we should be
    // prepared.
    CPPUNIT_ASSERT(MeshTools::valid_is_prepared(mesh));

    // Break some neighbor links.  Of course nobody would do this in
    // real life, right?
    Elem * elem0 = mesh.query_elem_ptr(0);
    if (elem0)
      for (auto n : elem0->side_index_range())
        {
          Elem * neigh = elem0->neighbor_ptr(n);
          if (!neigh)
            continue;
          auto n_neigh = neigh->which_neighbor_am_i(elem0);
          CPPUNIT_ASSERT(n_neigh != libMesh::invalid_uint);

          elem0->set_neighbor(n, nullptr);
          neigh->set_neighbor(n_neigh, nullptr);
        }

    // We're unprepared (prepare_for_use() will restitch those
    // neighbor pointers) but we're not marked that way.
    CPPUNIT_ASSERT(!MeshTools::valid_is_prepared(mesh));
  }

  void testDistributedMeshVerifyIsPrepared ()
  {
    DistributedMesh mesh(*TestCommWorld);
    testMeshBaseVerifyIsPrepared(mesh);
  }

  void testMeshVerifyIsPrepared ()
  {
    Mesh mesh(*TestCommWorld);
    testMeshBaseVerifyIsPrepared(mesh);
  }

  void testReplicatedMeshVerifyIsPrepared ()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    testMeshBaseVerifyIsPrepared(mesh);
  }
}; // End definition of class MeshBaseTest

CPPUNIT_TEST_SUITE_REGISTRATION( MeshBaseTest );
