// libMesh includes
#include <libmesh/libmesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>

// cppunit includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class AllSecondOrderTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllSecondOrderTest );

  CPPUNIT_TEST( allSecondOrder );
  CPPUNIT_TEST( allSecondOrderRange );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void allSecondOrder()
  {
    LOG_UNIT_TEST;

    DistributedMesh mesh(*TestCommWorld, /*dim=*/2);

    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_square(mesh, 2, 2);

    mesh.all_second_order();
  }

  void allSecondOrderRange()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    // Construct a single-element Quad4 mesh
    MeshTools::Generation::build_square(mesh, /*nx=*/1, /*ny=*/1);

    // Add extra nodes above the four original nodes
    Point z_trans(0,0,1);
    for (dof_id_type n=0; n<4; ++n)
      mesh.add_point(mesh.point(n) + z_trans, /*id=*/4+n, /*proc_id=*/0);

    // Construct Edge2 elements attached to Quad4 at individual points
    // and in subdomain 1.
    for (dof_id_type n=0; n<4; ++n)
      {
        Elem * elem = mesh.add_elem(new Edge2);
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(n+4);
        elem->subdomain_id() = 1;
      }

    // Convert only Edge2 elements (all elements in subdomain 1) to
    // SECOND-order.
    mesh.all_second_order_range(mesh.active_subdomain_elements_ptr_range(1),
                                /*full_ordered=*/true);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllSecondOrderTest );
