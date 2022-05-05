#include <libmesh/replicated_mesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/mesh_modification.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class CopyNodesAndElementsTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( CopyNodesAndElementsTest );

  CPPUNIT_TEST( replicatedCopy );
  CPPUNIT_TEST( distributedCopy );

  CPPUNIT_TEST_SUITE_END();

private:
  template <typename T>
  void collectMeshes()
  {
    // to trigger duplicate unique ids in the old code, you needed the following circumstances:
    // The single mesh generation routine (e.g. build_cube) *must* add a node after all the elements
    // such that it is not an element that has the largest unique id in the single mesh. (build_cube
    // with TET10 adds a node after all the elements are added becuase it calls all_second_order
    // after all the elements are created.) If an element has the highest unique id, then the final
    // value for _next_unique_id was (and still is) correct at the end of copy_nodes_and_elements
    // because we set the element unique ids before adding the elements into the original mesh, such
    // that _next_unique_id gets correctly updated. But when a node has the highest unique id, then
    // _next_unique_id would be incorrect with the old code (but is now correct with the new code!)
    // and after one mesh collection the value of parallel_max_unique_id() would be conceptually
    // wrong, such that the next time we collected a mesh we would be at risk of duplicating unique
    // ids
    T mesh(*TestCommWorld, /*dim=*/3);
    MeshTools::Generation::build_cube(mesh, 2, 2, 2, 0, 1, 0, 1, 0, 1, TET10);
    Real dmin = 2, dmax = 3;
    for (unsigned int num_other_meshes = 0; num_other_meshes < 2;
         ++num_other_meshes, dmin += 2, dmax += 2)
    {
      T other_mesh(*TestCommWorld, 3);
      MeshTools::Generation::build_cube(
          other_mesh, 2, 2, 2, dmin, dmax, dmin, dmax, dmin, dmax, TET10);
      dof_id_type node_delta = mesh.max_node_id();
      dof_id_type elem_delta = mesh.max_elem_id();
      unique_id_type unique_delta =
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          mesh.parallel_max_unique_id();
#else
          0;
#endif
      mesh.copy_nodes_and_elements(other_mesh, false, elem_delta, node_delta, unique_delta);
    }
  }

public:
  void setUp() {}

  void tearDown() {}

  void replicatedCopy() { LOG_UNIT_TEST; collectMeshes<ReplicatedMesh>(); }

  void distributedCopy() { LOG_UNIT_TEST; collectMeshes<DistributedMesh>(); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( CopyNodesAndElementsTest );
