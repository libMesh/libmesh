// libmesh includes
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh.h>

// unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class VolumeTest : public CppUnit::TestCase
{

public:
  CPPUNIT_TEST_SUITE( VolumeTest );
  CPPUNIT_TEST( testEdge3Volume );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testEdge3Volume()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line (mesh, /*nelem=*/1, 0., 1., EDGE3);
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(3), mesh.n_nodes());

    auto edge3 = mesh.elem_ptr(0);

    // Check unperturbed, straight edge case
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check middle node perturbed in +x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    Real & middle_node_x = edge3->node_ref(2)(0);
    middle_node_x += 1.e-3;
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check middle node perturbed in +x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    middle_node_x -= 2.e-3;
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( VolumeTest );
