#include <libmesh/libmesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class AllSecondOrderTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllSecondOrderTest );

  CPPUNIT_TEST( allSecondOrder );

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
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllSecondOrderTest );
