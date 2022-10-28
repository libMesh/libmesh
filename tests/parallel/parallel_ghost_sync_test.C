
// libMesh includes
#include <libmesh/parallel_ghost_sync.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


namespace {

  using namespace libMesh;

  struct SyncAndTest {
    typedef std::vector<Real> datum;

    const MeshBase & mesh;

    std::unordered_set<dof_id_type> answered_ids;

    SyncAndTest(MeshBase & _mesh) :
      mesh(_mesh) {}

    ~SyncAndTest() {
      for (const Node * node : mesh.node_ptr_range())
        if (node->processor_id() != mesh.processor_id())
          CPPUNIT_ASSERT(answered_ids.count(node->id()));
    }

    void gather_data (const std::vector<dof_id_type> & ids,
                      std::vector<datum> & data)
    {
      data.resize(ids.size());

      for (auto i : index_range(ids))
        {
          const Node * node = mesh.query_node_ptr(ids[i]);

          // We should only be asked about nodes we own
          CPPUNIT_ASSERT(node);
          CPPUNIT_ASSERT_EQUAL(mesh.processor_id(), node->processor_id());

          data[i] = {(*node)(2), (*node)(1), (*node)(0)};
        }
    }

    void act_on_data (const std::vector<dof_id_type> & ids,
                      const std::vector<datum> & data)
    {
      CPPUNIT_ASSERT_EQUAL(ids.size(), data.size());

      for (auto i : index_range(ids))
        {
          answered_ids.insert(ids[i]);

          // We should have every node we asked about
          const Node * node = mesh.query_node_ptr(ids[i]);
          CPPUNIT_ASSERT(node);
          CPPUNIT_ASSERT(mesh.processor_id() != node->processor_id());

          CPPUNIT_ASSERT_EQUAL(data[i].size(), std::size_t(3));
          CPPUNIT_ASSERT_EQUAL(data[i][0], (*node)(2));
          CPPUNIT_ASSERT_EQUAL(data[i][1], (*node)(1));
          CPPUNIT_ASSERT_EQUAL(data[i][2], (*node)(0));
        }
    }
  };
}


using namespace libMesh;

class ParallelGhostSyncTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ParallelGhostSyncTest );

  CPPUNIT_TEST( testByXYZ );

  CPPUNIT_TEST_SUITE_END();

private:
  std::unique_ptr<Mesh> _mesh;

public:
  void setUp()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);
    constexpr unsigned int n = 3;
    MeshTools::Generation::build_cube (*_mesh, n, n, n, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, HEX27);
  }

  void tearDown()
  {}

  void testByXYZ()
  {
    LOG_UNIT_TEST;

    SyncAndTest sync(*_mesh);
    LocationMap<Node> loc_map;
    loc_map.init(*_mesh);
    Parallel::sync_dofobject_data_by_xyz
      (*TestCommWorld, _mesh->nodes_begin(), _mesh->nodes_end(),
       loc_map, sync);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( ParallelGhostSyncTest );
