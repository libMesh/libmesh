// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/serial_mesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/remote_elem.h>

#include "test_comm.h"

using namespace libMesh;

class BoundaryMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh generates
   * boundary meshes correctly.
   */
public:
  CPPUNIT_TEST_SUITE( BoundaryMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

protected:

  Mesh* _mesh;
  Mesh* _all_boundary_mesh;
  Mesh* _left_boundary_mesh;

  void build_mesh()
  {
    _mesh = new Mesh(*TestCommWorld);
    _all_boundary_mesh = new Mesh(*TestCommWorld);
    _left_boundary_mesh = new Mesh(*TestCommWorld);

    MeshTools::Generation::build_square(*_mesh, 3, 5,
                                        0.1, 0.9, 0.1, 0.9, QUAD9);

    // We'll need to skip repartitioning with ParallelMesh for now;
    // otherwise the boundary meshes' interior parents might get
    // shuffled off to different processors.
    if (!_mesh->is_serial())
      {
        _mesh->skip_partitioning(true);
        _left_boundary_mesh->skip_partitioning(true);
        _all_boundary_mesh->skip_partitioning(true);
      }

    // Get the border of the square
    _mesh->boundary_info->sync(*_all_boundary_mesh);

    std::set<boundary_id_type> left_id, right_id;
    left_id.insert(3);
    right_id.insert(1);

    // Add the right side of the square to the square; this should
    // make it a mixed dimension mesh
    _mesh->boundary_info->add_elements(right_id, *_mesh);
    _mesh->prepare_for_use();

    // Add the left side of the square to its own boundary mesh.
    _mesh->boundary_info->sync(left_id, *_left_boundary_mesh);
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _all_boundary_mesh;
    delete _left_boundary_mesh;
    delete _mesh;
  }

  void testMesh()
  {
    // There'd better be 3*5 + 5 elements in the interior plus right
    // boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)20, _mesh->n_elem() );

    // There'd better be only 7*11 nodes in the interior
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)77, _mesh->n_nodes() );

    // There'd better be only 2*(3+5) elements on the full boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)16, _all_boundary_mesh->n_elem() );

    // There'd better be only 2*2*(3+5) nodes on the full boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)32, _all_boundary_mesh->n_nodes() );

    // There'd better be only 5 elements on the left boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)5, _left_boundary_mesh->n_elem() );

    // There'd better be only 2*5+1 nodes on the left boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)11, _left_boundary_mesh->n_nodes() );

    this->sanityCheck();
  }

  void sanityCheck()
  {
    // Sanity check all the elements
    MeshBase::const_element_iterator elem_it =
      _mesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end =
      _mesh->active_elements_end();
    for (; elem_it != elem_end; ++elem_it)
      {
        const Elem *elem = *elem_it;

        const Elem *pip = elem->interior_parent();

        // On a ParallelMesh we might not be able to see the
        // interior_parent of a non-local element
        if (pip == remote_elem)
          {
            CPPUNIT_ASSERT(elem->processor_id() != TestCommWorld->rank());
            continue;
          }

        // All the edges should have interior parents; none of the
        // quads should.
        if (pip)
          {
            CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);
            CPPUNIT_ASSERT_EQUAL(pip->type(), QUAD9);
            CPPUNIT_ASSERT_EQUAL(pip->level(), elem->level());

            // We only added right edges
            CPPUNIT_ASSERT_DOUBLES_EQUAL(elem->centroid()(0), 0.9,
                                         TOLERANCE*TOLERANCE);
          }
        else
          {
            CPPUNIT_ASSERT_EQUAL(elem->type(), QUAD9);
          }
      }

    MeshBase::const_element_iterator left_bdy_elem_it =
      _left_boundary_mesh->active_elements_begin();
    const MeshBase::const_element_iterator left_bdy_elem_end =
      _left_boundary_mesh->active_elements_end();
    for (; left_bdy_elem_it != left_bdy_elem_end; ++left_bdy_elem_it)
      {
        const Elem *elem = *left_bdy_elem_it;

        CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);

        const Elem *pip = elem->interior_parent();

        // On a ParallelMesh we might not be able to see the
        // interior_parent of a non-local element
        if (pip == remote_elem)
          {
            CPPUNIT_ASSERT(elem->processor_id() != TestCommWorld->rank());
            continue;
          }

        // All the edges should have interior parents
        CPPUNIT_ASSERT(pip);
        CPPUNIT_ASSERT_EQUAL(pip->type(), QUAD9);
        CPPUNIT_ASSERT_EQUAL(pip->level(), elem->level());

        // We only added left edges
        CPPUNIT_ASSERT_DOUBLES_EQUAL(elem->centroid()(0), 0.1,
                                     TOLERANCE*TOLERANCE);
      }


    MeshBase::const_element_iterator all_bdy_elem_it =
      _left_boundary_mesh->active_elements_begin();
    const MeshBase::const_element_iterator all_bdy_elem_end =
      _left_boundary_mesh->active_elements_end();
    for (; all_bdy_elem_it != all_bdy_elem_end; ++all_bdy_elem_it)
      {
        const Elem *elem = *all_bdy_elem_it;

        CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);

        const Elem *pip = elem->interior_parent();

        // On a ParallelMesh we might not be able to see the
        // interior_parent of a non-local element
        if (pip == remote_elem)
          {
            CPPUNIT_ASSERT(elem->processor_id() != TestCommWorld->rank());
            continue;
          }

        // All the edges should have interior parents
        CPPUNIT_ASSERT(pip);
        CPPUNIT_ASSERT_EQUAL(pip->type(), QUAD9);
        CPPUNIT_ASSERT_EQUAL(pip->level(), elem->level());
      }

  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryMeshTest );

#ifdef LIBMESH_ENABLE_AMR

class BoundaryRefinedMeshTest : public BoundaryMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  CPPUNIT_TEST_SUITE( BoundaryRefinedMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
    this->build_mesh();

    // Need to refine interior mesh before separate boundary meshes,
    // if we want to get interior_parent links right.
    MeshRefinement(*_mesh).uniformly_refine(1);
    MeshRefinement(*_left_boundary_mesh).uniformly_refine(1);
    MeshRefinement(*_all_boundary_mesh).uniformly_refine(1);
  }

  void testMesh()
  {
    // There'd better be 3*5*4 + 5*2 active elements in the interior
    // plus right boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)70, _mesh->n_active_elem() );

    // Plus the original 20 now-inactive elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)90, _mesh->n_elem() );

    // There'd better be only 13*21 nodes in the interior
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)273, _mesh->n_nodes() );

    // There'd better be only 2*2*(3+5) active elements on the full boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)32,
                          _all_boundary_mesh->n_active_elem() );

    // Plus the original 16 now-inactive elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)48, _all_boundary_mesh->n_elem() );

    // There'd better be only 2*2*2*(3+5) nodes on the full boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)64, _all_boundary_mesh->n_nodes() );

    // There'd better be only 2*5 active elements on the left boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)10,
                          _left_boundary_mesh->n_active_elem() );

    // Plus the original 5 now-inactive elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)15, _left_boundary_mesh->n_elem() );

    // There'd better be only 2*2*5+1 nodes on the left boundary
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)21, _left_boundary_mesh->n_nodes() );

    this->sanityCheck();
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryRefinedMeshTest );

#endif
