#include "libmesh/cell_c0polyhedron.h"
#include "libmesh/checkpoint_io.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/face_c0polygon.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/replicated_mesh.h"

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class CheckpointIOTest : public CppUnit::TestCase {
  /**
   * This test verifies that we can write files with the CheckpointIO object.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( CheckpointIOTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testAsciiDistRepSplitter );
  CPPUNIT_TEST( testBinaryDistRepSplitter );
  CPPUNIT_TEST( testAsciiRepDistSplitter );
  CPPUNIT_TEST( testBinaryRepDistSplitter );
  CPPUNIT_TEST( testAsciiRepRepSplitter );
  CPPUNIT_TEST( testBinaryRepRepSplitter );
  CPPUNIT_TEST( testAsciiDistDistSplitter );
  CPPUNIT_TEST( testBinaryDistDistSplitter );
  CPPUNIT_TEST( testAsciiC0Polygon );
  CPPUNIT_TEST( testBinaryC0Polygon );
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testAsciiC0Polyhedron );
  CPPUNIT_TEST( testBinaryC0Polyhedron );
#endif
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
  void testSplitter(bool binary, bool using_distmesh, bool skip_partition = false)
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
      if (skip_partition)
        mesh.skip_partitioning(true);

      CheckpointIO cpr(mesh);
      cpr.current_n_processors() = n_procs;
      cpr.binary() = binary;
      cpr.read(filename);

      // If we decided to skip partitioning, then we shouldn't be
      // waiting for a partition() call to get our n_partitions()
      // cache in order
      if (skip_partition)
        CPPUNIT_ASSERT_EQUAL(mesh.n_partitions(), n_procs);

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

  void testC0PolygonCheckpoint(bool binary)
  {
#ifdef LIBMESH_HAVE_XDR
    // Parallel element packing does not yet support runtime topology,
    // so exercise serial CheckpointIO independently on rank 0.
    if (TestCommWorld->rank() != 0)
      return;

    Parallel::Communicator comm_self;
    const std::string filename =
      std::string("checkpoint_c0polygon.cp") + (binary ? "r" : "a");
    const std::vector<Point> points =
      {{0., 0.}, {1., 0.}, {1.5, 0.5}, {1., 1.}, {0., 1.}};

    {
      Mesh mesh(comm_self, 2);
      mesh.allow_renumbering(false);

      auto polygon =
        std::make_unique<C0Polygon>(cast_int<unsigned int>(points.size()));
      polygon->set_id() = 0;
      for (auto n : index_range(points))
        polygon->set_node(n, mesh.add_point(points[n], n));

      mesh.add_elem(std::move(polygon));
      mesh.prepare_for_use();

      CheckpointIO checkpoint(mesh, binary);
      checkpoint.write(filename);
    }

    {
      Mesh mesh(comm_self);
      CheckpointIO checkpoint(mesh, binary);
      checkpoint.read(filename);

      CPPUNIT_ASSERT_EQUAL(dof_id_type(1), mesh.n_elem());
      CPPUNIT_ASSERT_EQUAL(dof_id_type(5), mesh.n_nodes());
      CPPUNIT_ASSERT_EQUAL(2u, static_cast<unsigned int>(mesh.mesh_dimension()));

      const Elem * elem = mesh.query_elem_ptr(0);
      bool found_elem = elem;
      mesh.comm().max(found_elem);
      CPPUNIT_ASSERT(found_elem);

      if (elem)
        {
          CPPUNIT_ASSERT_EQUAL(C0POLYGON, elem->type());
          CPPUNIT_ASSERT_EQUAL(5u, elem->n_nodes());
          CPPUNIT_ASSERT_EQUAL(5u, elem->n_vertices());
          CPPUNIT_ASSERT_EQUAL(5u, elem->n_sides());
          CPPUNIT_ASSERT_EQUAL(5u, elem->n_edges());

          for (auto n : index_range(points))
            {
              CPPUNIT_ASSERT_EQUAL(cast_int<dof_id_type>(n), elem->node_id(n));
              CPPUNIT_ASSERT_EQUAL(points[n], elem->point(n));

              const auto side_nodes = elem->nodes_on_side(n);
              CPPUNIT_ASSERT_EQUAL(std::size_t(2), side_nodes.size());
              CPPUNIT_ASSERT_EQUAL(cast_int<unsigned int>(n), side_nodes[0]);
              CPPUNIT_ASSERT_EQUAL
                (cast_int<unsigned int>((n + 1) % points.size()), side_nodes[1]);
            }
        }
    }

    CheckpointIO::cleanup(filename, 1);
#endif // LIBMESH_HAVE_XDR
  }

#if LIBMESH_DIM > 2
  void testC0PolyhedronCheckpoint(bool binary)
  {
#ifdef LIBMESH_HAVE_XDR
    // Parallel element packing does not yet support runtime topology,
    // so exercise serial CheckpointIO independently on rank 0.
    if (TestCommWorld->rank() != 0)
      return;

    Parallel::Communicator comm_self;
    const std::string filename =
      std::string("checkpoint_c0polyhedron.cp") + (binary ? "r" : "a");
    const std::vector<Point> points =
      {{0., -2., 0.}, {-1., -1., 0.}, {-1., 1., 0.},
       {0., 2., 0.}, {1., 1., 0.}, {1., -1., 0.},
       {0., -2., 1.}, {-1., -1., 1.}, {-1., 1., 1.},
       {0., 2., 1.}, {1., 1., 1.}, {1., -1., 1.}};
    const std::vector<std::vector<unsigned int>> nodes_on_sides =
      {{0, 1, 2, 3, 4, 5},
       {0, 1, 7, 6},
       {1, 2, 8, 7},
       {2, 3, 9, 8},
       {3, 4, 10, 9},
       {4, 5, 11, 10},
       {5, 0, 6, 11},
       {6, 7, 8, 9, 10, 11}};

    {
      Mesh mesh(comm_self, 3);
      mesh.allow_renumbering(false);
      for (auto n : index_range(points))
        mesh.add_point(points[n], n);

      std::vector<std::shared_ptr<Polygon>> sides(nodes_on_sides.size());
      for (auto s : index_range(nodes_on_sides))
        {
          const auto & side_nodes = nodes_on_sides[s];
          auto side =
            std::make_shared<C0Polygon>
              (cast_int<unsigned int>(side_nodes.size()));
          for (auto n : index_range(side_nodes))
            side->set_node(n, mesh.node_ptr(side_nodes[n]));
          sides[s] = std::move(side);
        }

      std::unique_ptr<Node> mid_elem_node;
      auto polyhedron =
        std::make_unique<C0Polyhedron>(sides, mid_elem_node);
      CPPUNIT_ASSERT(mid_elem_node);

      // Explicit id so DistributedMesh matches ReplicatedMesh
      mid_elem_node->set_id(points.size());

      mesh.add_node(std::move(mid_elem_node));
      polyhedron->set_id() = 0;
      mesh.add_elem(std::move(polyhedron));
      mesh.prepare_for_use();

      CheckpointIO checkpoint(mesh, binary);
      checkpoint.write(filename);
    }

    {
      Mesh mesh(comm_self);
      CheckpointIO checkpoint(mesh, binary);
      checkpoint.read(filename);

      CPPUNIT_ASSERT_EQUAL(dof_id_type(1), mesh.n_elem());
      CPPUNIT_ASSERT_EQUAL(dof_id_type(13), mesh.n_nodes());
      CPPUNIT_ASSERT_EQUAL(3u, static_cast<unsigned int>(mesh.mesh_dimension()));

      const Elem * elem = mesh.query_elem_ptr(0);
      bool found_elem = elem;
      mesh.comm().max(found_elem);
      CPPUNIT_ASSERT(found_elem);

      if (elem)
        {
          CPPUNIT_ASSERT_EQUAL(C0POLYHEDRON, elem->type());
          CPPUNIT_ASSERT_EQUAL(12u, elem->n_vertices());
          CPPUNIT_ASSERT_EQUAL(13u, elem->n_nodes());
          CPPUNIT_ASSERT_EQUAL(8u, elem->n_sides());
          CPPUNIT_ASSERT_EQUAL(18u, elem->n_edges());
          CPPUNIT_ASSERT_EQUAL(dof_id_type(12),
                               elem->node_id(elem->n_vertices()));

          for (auto s : index_range(nodes_on_sides))
            {
              const auto side_nodes = elem->nodes_on_side(s);
              CPPUNIT_ASSERT_EQUAL(nodes_on_sides[s].size(), side_nodes.size());
              for (auto n : index_range(side_nodes))
                CPPUNIT_ASSERT_EQUAL
                  (cast_int<dof_id_type>(nodes_on_sides[s][n]),
                   elem->node_id(side_nodes[n]));
            }
        }
    }

    CheckpointIO::cleanup(filename, 1);
#endif // LIBMESH_HAVE_XDR
  }
#endif

  void testAsciiDistRepSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<DistributedMesh, ReplicatedMesh>(false, true);
  }

  void testBinaryDistRepSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<DistributedMesh, ReplicatedMesh>(true, true);
  }

  void testAsciiRepDistSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<ReplicatedMesh, DistributedMesh>(false, true);
  }

  void testBinaryRepDistSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<ReplicatedMesh, DistributedMesh>(true, true);
  }

  void testAsciiRepRepSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<ReplicatedMesh, ReplicatedMesh>(false, false);
  }

  void testBinaryRepRepSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<ReplicatedMesh, ReplicatedMesh>(true, false);
  }

  void testAsciiDistDistSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<DistributedMesh, DistributedMesh>(false, true);
  }

  void testBinaryDistDistSplitter()
  {
    LOG_UNIT_TEST;

    testSplitter<DistributedMesh, DistributedMesh>(true, true);
  }

  void testAsciiC0Polygon()
  {
    LOG_UNIT_TEST;

    testC0PolygonCheckpoint(false);
  }

  void testBinaryC0Polygon()
  {
    LOG_UNIT_TEST;

    testC0PolygonCheckpoint(true);
  }

#if LIBMESH_DIM > 2
  void testAsciiC0Polyhedron()
  {
    LOG_UNIT_TEST;

    testC0PolyhedronCheckpoint(false);
  }

  void testBinaryC0Polyhedron()
  {
    LOG_UNIT_TEST;

    testC0PolyhedronCheckpoint(true);
  }
#endif

  void testAsciiDistDistSplitterCache()
  {
    LOG_UNIT_TEST;

    testSplitter<DistributedMesh, DistributedMesh>(false, true, true);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( CheckpointIOTest );
