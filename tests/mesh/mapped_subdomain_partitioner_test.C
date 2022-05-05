#include <libmesh/libmesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include "libmesh/mapped_subdomain_partitioner.h"

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MappedSubdomainPartitionerTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * MappedSubdomainPartitioner on different numbers of processors.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MappedSubdomainPartitionerTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMappedSubdomainPartitioner );
#endif

  CPPUNIT_TEST_SUITE_END();

  /**
   * Note: this second public is necessary, something in the macros
   * above leaves us in a private region.
   */
public:
  void setUp() {}

  void tearDown() {}

  void testMappedSubdomainPartitioner()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    Real
      xmin = 0., xmax = 1.,
      ymin = 0., ymax = 10.;

    MeshTools::Generation::build_square (mesh,
                                         /*nx*/10,
                                         /*ny*/100,
                                         xmin, xmax,
                                         ymin, ymax,
                                         QUAD4);

    // The MappedSubdomainPartitioner partitions based on user-defined
    // assignment of subdomains to processors.
    mesh.partitioner() = std::make_unique<MappedSubdomainPartitioner>();

    // Get a pointer to the MappedSubdomainPartitioner so we can call its
    // API specifically.
    MappedSubdomainPartitioner * subdomain_partitioner =
      dynamic_cast<MappedSubdomainPartitioner *>(mesh.partitioner().get());

    // Create 2x as many subdomains as processors, then assign them in
    // the following way:
    // subdomains(0,1) -> processor 0
    // subdomains(2,3) -> processor 1
    // subdomains(4,5) -> processor 2
    // ...
    // subdomains(n,n+1) -> processor n/2
    subdomain_id_type n_subdomains = 2 * TestCommWorld->size();
    for (subdomain_id_type sbd_id=0; sbd_id<n_subdomains; sbd_id+=2)
      {
        subdomain_partitioner->subdomain_to_proc[sbd_id] = sbd_id/2;
        subdomain_partitioner->subdomain_to_proc[sbd_id+1] = sbd_id/2;
      }

    // Assign subdomain ids to elements sequentially.
    {
      subdomain_id_type current_subdomain_id = 0;
      for (auto & elem : mesh.element_ptr_range())
        {
          elem->subdomain_id() = current_subdomain_id++;

          // Wrap around
          if (current_subdomain_id == n_subdomains)
            current_subdomain_id = 0;
        }
    }

    // Partition again, now that we have set up the MappedSubdomainPartitioner.
    mesh.partition();

    // Assert that the partitioning worked as expected.
    for (auto & elem : mesh.element_ptr_range())
      {
        // Subdomain id n should map to processor id n/2.
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(elem->subdomain_id()/2),
                             static_cast<int>(elem->processor_id()));
      }
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MappedSubdomainPartitionerTest );
