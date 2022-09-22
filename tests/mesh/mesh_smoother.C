#include <libmesh/replicated_mesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_smoother_laplace.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/exodusII_io.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class SmoothMeshTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(SmoothMeshTest);

  CPPUNIT_TEST(laplace_smooth_2d_fixed_boundary_nodes);
  CPPUNIT_TEST(laplace_smooth_2d_movable_boundary_nodes);
  CPPUNIT_TEST(laplace_smooth_3d_fixed_boundary_nodes);
  CPPUNIT_TEST(laplace_smooth_3d_movable_boundary_nodes);

  CPPUNIT_TEST_SUITE_END();

private:
  template <typename T>
  void test(const std::string in_file, const std::string gold_file, bool move_nodes)
  {
    // load input mesh
    T in_mesh(*TestCommWorld);
    ExodusII_IO in_exii(in_mesh);
    if (in_mesh.processor_id() == 0)
      in_exii.read(in_file);
    MeshCommunication().broadcast(in_mesh);
    in_mesh.skip_partitioning(true);
    in_mesh.allow_renumbering(false);
    in_mesh.prepare_for_use();

    // load gold file
    T gold_mesh(*TestCommWorld);
    ExodusII_IO gold_exii(gold_mesh);
    if (gold_mesh.processor_id() == 0)
      gold_exii.read(gold_file);
    MeshCommunication().broadcast(gold_mesh);
    gold_mesh.skip_partitioning(true);
    gold_mesh.allow_renumbering(false);
    gold_mesh.prepare_for_use();

    // run mesh smoother for 100 iterations
    LaplaceMeshSmoother lms(static_cast<UnstructuredMesh &>(in_mesh), move_nodes);
    lms.smooth(100);

    // compare node positions
    for (const auto &node : in_mesh.node_ptr_range())
    {
      Point a = *node;
      Point b = gold_mesh.node_ref(node->id());
      // std::cout << node->id() << ": " << a << " == " << b << " (" << (a - b).norm() << ")\n";
      CPPUNIT_ASSERT((a - b).norm() < 1e-6);
    }
  }

public:
  void setUp() {}

  void tearDown() {}

  void laplace_smooth_2d_fixed_boundary_nodes()
  {
    LOG_UNIT_TEST;
    test<DistributedMesh>("meshes/smooth2d_in.e", "meshes/smooth2d_nomove.e", false);
    test<ReplicatedMesh>("meshes/smooth2d_in.e", "meshes/smooth2d_nomove.e", false);
  }
  void laplace_smooth_2d_movable_boundary_nodes()
  {
    LOG_UNIT_TEST;
    test<DistributedMesh>("meshes/smooth2d_in.e", "meshes/smooth2d_move.e", true);
    test<ReplicatedMesh>("meshes/smooth2d_in.e", "meshes/smooth2d_move.e", true);
  }
  void laplace_smooth_3d_fixed_boundary_nodes()
  {
    LOG_UNIT_TEST;
    test<DistributedMesh>("meshes/smooth3d_in.e", "meshes/smooth3d_nomove.e", false);
    test<ReplicatedMesh>("meshes/smooth3d_in.e", "meshes/smooth3d_nomove.e", false);
  }
  void laplace_smooth_3d_movable_boundary_nodes()
  {
    LOG_UNIT_TEST;
    test<DistributedMesh>("meshes/smooth3d_in.e", "meshes/smooth3d_move.e", true);
    test<ReplicatedMesh>("meshes/smooth3d_in.e", "meshes/smooth3d_move.e", true);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(SmoothMeshTest);
