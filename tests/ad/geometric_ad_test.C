#include <libmesh/distributed_mesh.h>
#include <libmesh/instantiate_real_type.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_METAPHYSICL

#include "metaphysicl/dualdynamicsparsenumberarray.h"

typedef MetaPhysicL::DualNumber<Real, MetaPhysicL::DynamicSparseNumberArray<Real, unsigned int>>
DualReal;

namespace libMesh
{
INSTANTIATE_ALL_REAL_TYPE(DualReal);
}

using namespace libMesh;


class GeometricADTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(GeometricADTest);

  CPPUNIT_TEST(testMechanicalContact);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testMechanicalContact()
  {

    // This is the reference mesh
    DistributedMesh mesh(*TestCommWorld);

    // This is the displaced mesh
    DistributedMeshTempl<DualReal> disp_mesh(*TestCommWorld);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(GeometricADTest);

#endif // LIBMESH_HAVE_METAPHYSICL
