#ifndef __fe_test_h__
#define __fe_test_h__

#include "test_comm.h"

#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/partitioner.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>

// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

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

    // Unfortunately, it turns out that our METIS and ParMETIS
    // partitioners *do* suck, and can't reliabily give us more than
    // 13 non-empty ranks on the above 27 element mesh.
    CPPUNIT_ASSERT(nonempty_procs >= n_nonempty ||
                   nonempty_procs >= 13);
  }

  void testPartitionEmpty()
  {
    MeshClass mesh(*TestCommWorld);
    PartitionerSubclass newpart;

    // With a 0 element mesh this should just give us 0 subpartitions
    // regardless of n_procs
    newpart.partition(mesh, TestCommWorld->size());
  }

  void testPartition1()
  {
    this->testPartition(1);
  }

  void testPartition2()
  {
    this->testPartition(2);
  }

  void testPartitionNProc()
  {
    this->testPartition(TestCommWorld->size());
  }
};

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We'll put an
// ignore_warnings at the end of this file so it's the last warnings
// related header that our including code sees.
#include <libmesh/ignore_warnings.h>

#define INSTANTIATE_PARTITIONER_TEST(partitionersubclass, meshclass)        \
  class PartitionerTest_##partitionersubclass##_##meshclass :               \
    public PartitionerTest<partitionersubclass, meshclass> {                \
  public:                                                                   \
  CPPUNIT_TEST_SUITE( PartitionerTest_##partitionersubclass##_##meshclass); \
  PARTITIONERTEST                                                           \
  CPPUNIT_TEST_SUITE_END();                                                 \
  };                                                                        \
                                                                            \
  CPPUNIT_TEST_SUITE_REGISTRATION( PartitionerTest_##partitionersubclass##_##meshclass );

#endif // #ifdef __fe_test_h__
