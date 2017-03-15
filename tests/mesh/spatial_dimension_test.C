// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

class MeshSpatialDimensionTest : public CppUnit::TestCase
{
  /**
   * This unit test tests a few basic things with the new
   * Mesh::spatial_dimension() implementation.
   */
public:
  CPPUNIT_TEST_SUITE( MeshSpatialDimensionTest );

  CPPUNIT_TEST( test1D );
  CPPUNIT_TEST( test2D );

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
    // 1.) Test that build_line() produces a Mesh with spatial_dimension==1
    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line (mesh, /*n_elem=*/2, /*xmin=*/0., /*xmax=*/1., EDGE2);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.spatial_dimension());

    // 2.) Move the nodes in the y-direction, test that spatial_dimension==2
    // The spatial dimension is updated during prepare_for_use().
    {
      MeshBase::node_iterator node_it  = mesh.nodes_begin();
      MeshBase::node_iterator node_end = mesh.nodes_end();

      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;
          node(1) = node(0)*node(0);
        }
    }

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());


    // 3.) Move nodes back to zero, check that spatial_dimension is *not* decreased
    {
      MeshBase::node_iterator node_it  = mesh.nodes_begin();
      MeshBase::node_iterator node_end = mesh.nodes_end();

      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;
          node(1) = 0.;
        }
    }
    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.spatial_dimension());


    // 4.) Move z-coordinate of nodes, check that spatial_dimension is now 3.
#if LIBMESH_DIM == 3
    {
      MeshBase::node_iterator node_it  = mesh.nodes_begin();
      MeshBase::node_iterator node_end = mesh.nodes_end();

      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;
          node(2) = node(0)*node(0);
        }
    }

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), mesh.spatial_dimension());
#endif
  }



  void test2D()
  {
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
    {
      MeshBase::node_iterator node_it  = mesh.nodes_begin();
      MeshBase::node_iterator node_end = mesh.nodes_end();

      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;
          node(2) = node(0)*node(0) + node(1)*node(1);
        }
    }

    mesh.prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), mesh.mesh_dimension());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(3), mesh.spatial_dimension());


    // 3.) Move nodes back to zero, check that spatial_dimension is *not* decreased
    {
      MeshBase::node_iterator node_it  = mesh.nodes_begin();
      MeshBase::node_iterator node_end = mesh.nodes_end();

      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;
          node(2) = 0.;
        }
    }
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
