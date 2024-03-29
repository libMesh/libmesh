#include <libmesh/equation_systems.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshSpatialDimensionTest : public CppUnit::TestCase
{
  /**
   * This unit test tests a few basic things with the new
   * Mesh::spatial_dimension() implementation.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshSpatialDimensionTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( test1D );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( test2D );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void test1D()
  {
    LOG_UNIT_TEST;

    // 1.) Test that build_line() produces a Mesh with spatial_dimension==1
    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line (mesh, /*n_elem=*/2, /*xmin=*/0., /*xmax=*/1., EDGE2);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.spatial_dimension());

    // 2.) Move the nodes in the y-direction, test that spatial_dimension==2
    // The spatial dimension is updated during prepare_for_use().
    for (auto & node : mesh.node_ptr_range())
      (*node)(1) = (*node)(0) * (*node)(0);

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());


    // 3.) Move nodes back to zero, check that spatial_dimension is *not* decreased
    for (auto & node : mesh.node_ptr_range())
      (*node)(1) = 0.;

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());


    // 4.) Move z-coordinate of nodes, check that spatial_dimension is now 3.
#if LIBMESH_DIM == 3
    for (auto & node : mesh.node_ptr_range())
      (*node)(2) = (*node)(0) * (*node)(0);

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), mesh.spatial_dimension());
#endif
  }



  void test2D()
  {
    LOG_UNIT_TEST;

    // 1.) Test that build_cube() produces a Mesh with spatial_dimension==2
    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square (mesh,
                                         /*nx=*/2, /*ny=*/2,
                                         /*xmin=*/0., /*xmax=*/1.,
                                         /*ymin=*/0., /*ymax=*/1.,
                                         QUAD4);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());

    // 2.) Move the nodes in the z-direction, test that spatial_dimension==3
    // The spatial dimension is updated during prepare_for_use().
    for (auto & node : mesh.node_ptr_range())
      (*node)(2) =
        (*node)(0) * (*node)(0) +
        (*node)(1) * (*node)(1);

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), mesh.spatial_dimension());


    // 3.) Move nodes back to zero, check that spatial_dimension is *not* decreased
    for (auto & node : mesh.node_ptr_range())
      (*node)(2) = 0.;

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), mesh.spatial_dimension());

    // 4.) The user is allowed to set the spatial dimension, make sure
    // that's honored if it's consistent.
    mesh.set_spatial_dimension(2);
    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());

    // 5.) If the user sets a spatial dimension that's not big enough, make sure that
    // prepare_for_use() increases it.
    mesh.set_spatial_dimension(1);
    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshSpatialDimensionTest );
