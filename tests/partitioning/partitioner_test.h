#ifndef __fe_test_h__
#define __fe_test_h__

#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/partitioner.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


#if LIBMESH_DIM > 2
#define PARTITIONERTEST                         \
  CPPUNIT_TEST( testPartitionEmpty );           \
  CPPUNIT_TEST( testPartition1 );               \
  CPPUNIT_TEST( testPartition2 );               \
  CPPUNIT_TEST( testPartitionNProc );
#else
#define PARTITIONERTEST
#endif

using namespace libMesh;

template <typename PartitionerSubclass, typename MeshClass>
class PartitionerTest : public CppUnit::TestCase {
protected:

  std::string libmesh_suite_name;

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testPartition(processor_id_type n_parts)
  {
    MeshClass mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1., 0., 1., 0., 1.,
                                       HEX8);

    // FIXME: Splitting meshes into more than n_proc parts currently
    // requires us to start with a mesh entirely assigned to proc 0
    PartitionerSubclass newpart;
    newpart.partition(mesh, 1);

    for (auto elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(elem->processor_id(), processor_id_type(0));

    for (auto node : mesh.node_ptr_range())
      CPPUNIT_ASSERT_EQUAL(node->processor_id(), processor_id_type(0));

    // But then we can manually partition, into at most many parts as
    // were requested.
    newpart.partition(mesh, n_parts);

    // We expect the partitioner not to suck - every processor
    // rank (up to n_elem()) ought to have at least one element on it.
    const processor_id_type n_nonempty =
      std::min(n_parts, processor_id_type(3*3*3));

    // Let's make sure we can see them all even on a DistributedMesh
    mesh.allgather();

    processor_id_type nonempty_procs = 0;
    for (processor_id_type p=0; p != n_nonempty; ++p)
      {
        const std::size_t n_elem_on_p =
          std::distance(mesh.pid_elements_begin(p),
                        mesh.pid_elements_end(p));
        if (n_elem_on_p)
          nonempty_procs++;
      }

    // Unfortunately, it turns out that our partitioners *do* suck:
    //
    // * Metis and ParMetis can't reliabily give us more than
    // 13 non-empty ranks on the above 27 element mesh.
    //
    // * Our SFC partitioners are suboptimal with 8 or more ranks for
    // only 27 elements: the smallest integer block size for 27/8
    // elements-per-rank that won't overflow is 4, but that only fills
    // 7 ranks, because we don't equidistribute the underflow.
    CPPUNIT_ASSERT(nonempty_procs >= n_nonempty ||
                   nonempty_procs >= 13 ||
                   (nonempty_procs >= 7 &&
                    n_nonempty == 8) ||
                   (nonempty_procs >= 9 &&
                    n_nonempty >= 10));
  }

  void testPartitionEmpty()
  {
    LOG_UNIT_TEST;

    MeshClass mesh(*TestCommWorld);
    PartitionerSubclass newpart;

    // With a 0 element mesh this should just give us 0 subpartitions
    // regardless of n_procs
    newpart.partition(mesh, TestCommWorld->size());
  }

  void testPartition1()
  {
    LOG_UNIT_TEST;

    this->testPartition(1);
  }

  void testPartition2()
  {
    LOG_UNIT_TEST;

    this->testPartition(2);
  }

  void testPartitionNProc()
  {
    LOG_UNIT_TEST;

    this->testPartition(TestCommWorld->size());
  }
};

#define INSTANTIATE_PARTITIONER_TEST(partitionersubclass, meshclass)        \
  class PartitionerTest_##partitionersubclass##_##meshclass :               \
    public PartitionerTest<partitionersubclass, meshclass> {                \
  public:                                                                   \
  PartitionerTest_##partitionersubclass##_##meshclass() :                   \
    PartitionerTest<partitionersubclass,meshclass>() {                      \
    if (unitlog->summarized_logs_enabled())                                 \
      this->libmesh_suite_name = "PartitionerTest";                         \
    else                                                                    \
      this->libmesh_suite_name = "PartitionerTest_" #partitionersubclass "_" #meshclass; \
  }                                                                         \
  CPPUNIT_TEST_SUITE( PartitionerTest_##partitionersubclass##_##meshclass); \
  PARTITIONERTEST                                                           \
  CPPUNIT_TEST_SUITE_END();                                                 \
  };                                                                        \
                                                                            \
  CPPUNIT_TEST_SUITE_REGISTRATION( PartitionerTest_##partitionersubclass##_##meshclass );

#endif // #ifdef __fe_test_h__
