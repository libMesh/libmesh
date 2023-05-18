#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include "libmesh/enum_partitioner_type.h"
#include <libmesh/mesh_generation.h>
#include <libmesh/parallel.h>
#include <libmesh/partitioner.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class MeshDeletionsTest : public CppUnit::TestCase {
  /**
   * This test ensures operations like delete_elem don't cause
   * problems with later mesh preparation
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshDeletionsTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testDeleteElemReplicated );
  CPPUNIT_TEST( testDeleteElemDistributed );
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

  void testDeleteElem(UnstructuredMesh & mesh)
  {
    LOG_UNIT_TEST;

    // Let's not trust Metis or Parmetis to give us a replicable test
    // case.
    mesh.partitioner() = Partitioner::build(LINEAR_PARTITIONER);

    MeshTools::Generation::build_square(mesh,
                                        3, 3,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    std::vector<std::vector<int>> deletions =
      {{1, 2, 4}, {0, 3}, {0, 2}, {1}};
    dof_id_type n_elem = 9;
    for (auto deletionset : deletions)
      {
        for (int e : deletionset)
          {
            --n_elem;
            Elem * elem = mesh.query_elem_ptr(e);
            if (elem)
              {
                elem->remove_links_to_me();
                mesh.delete_elem(elem);
              }
          }

        mesh.prepare_for_use();

        CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_elem);
      }
  }

  void testDeleteElemReplicated()
  {
    LOG_UNIT_TEST;
    ReplicatedMesh mesh(*TestCommWorld);
    testDeleteElem(mesh);
  }

  void testDeleteElemDistributed()
  {
    LOG_UNIT_TEST;
    DistributedMesh mesh(*TestCommWorld);
    testDeleteElem(mesh);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshDeletionsTest );
