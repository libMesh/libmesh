#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_serializer.h>
#include <libmesh/remote_elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <memory>

using namespace libMesh;

class BoundaryMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh generates
   * boundary meshes correctly.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( BoundaryMeshTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testBoundarySerialization );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  std::unique_ptr<UnstructuredMesh> _mesh;
  std::unique_ptr<UnstructuredMesh> _multi_boundary_mesh;
  std::unique_ptr<UnstructuredMesh> _exterior_boundary_mesh;
  std::unique_ptr<UnstructuredMesh> _left_boundary_mesh;
  std::unique_ptr<UnstructuredMesh> _internal_boundary_mesh;

  // We create internal sideset IDs that do not conflict with
  // sidesets 0-3 that get created by build_square().
  static constexpr boundary_id_type bid1 = 5;
  static constexpr boundary_id_type bid2 = 6;

  void build_mesh()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);
    _multi_boundary_mesh = std::make_unique<Mesh>(*TestCommWorld);
    _exterior_boundary_mesh = std::make_unique<Mesh>(*TestCommWorld);

    // We want to test Distributed->Replicated sync; this does that in
    // some builds
    _left_boundary_mesh = std::make_unique<ReplicatedMesh>(*TestCommWorld);

    // We want to test Replicated->Distributed sync; this does that in
    // other builds
    _internal_boundary_mesh = std::make_unique<DistributedMesh>(*TestCommWorld);

    MeshTools::Generation::build_square(*_mesh, 3, 5,
                                        0.2, 0.8, 0.2, 0.7, QUAD9);

    // We'll need to skip most repartitioning with DistributedMesh for
    // now; otherwise the boundary meshes' interior parents might get
    // shuffled off to different processors.
    if (!_mesh->is_serial())
      {
        _mesh->skip_noncritical_partitioning(true);
        _left_boundary_mesh->skip_noncritical_partitioning(true);
        _multi_boundary_mesh->skip_noncritical_partitioning(true);
        _exterior_boundary_mesh->skip_noncritical_partitioning(true);
        _internal_boundary_mesh->skip_noncritical_partitioning(true);
      }

    // Set subdomain ids for specific elements. This allows us to later
    // build an internal sideset with respect to a given
    // subdomain. The element subdomains look like:
    // ___________________
    // |  2  |  2  |  2  |
    // |_____|_____|_____|
    // |  2  |  2  |  2  |
    // |_____|_____|_____|
    // |  2  |  2  |  2  |
    // |_____|_____|_____|
    // |  1  |  1  |  2  |
    // |_____|_____|_____|
    // |  1  |  1  |  2  |
    // |_____|_____|_____|
    //
    // and we will create an internal sideset along the border between
    // subdomains 1 and 2.

    for (auto & elem : _mesh->active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        if (c(0) < 0.6 && c(1) < 0.4)
          elem->subdomain_id() = 1;
        else
          elem->subdomain_id() = 2;
      }

    // To test the "relative to" feature, we add the same sides to the
    // same sideset twice, from elements in subdomain 2 the second
    // time.  These should not show up in the BoundaryMesh, i.e. there
    // should not be overlapped elems in the BoundaryMesh.
    BoundaryInfo & bi = _mesh->get_boundary_info();

    for (auto & elem : _mesh->active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        if (c(0) < 0.6 && c(1) < 0.4)
          {
            if (c(0) > 0.4)
              bi.add_side(elem, 1, bid1);
            if (c(1) > 0.3)
              bi.add_side(elem, 2, bid1);
          }
        else
          {
            if (c(0) < 0.75 && c(1) < 0.4)
              bi.add_side(elem, 3, bid2);
            if (c(0) < 0.6 && c(1) < 0.5)
              bi.add_side(elem, 0, bid2);
          }
      }
  }

  void sync_and_test_meshes()
  {
    // Get the border of the square
    const std::set<boundary_id_type> exterior_boundaries {0, 1, 2, 3};
    _mesh->get_boundary_info().sync(exterior_boundaries,
                                    *_exterior_boundary_mesh);

    // The mesh of all boundaries is, because of the two nodes each
    // joining three edges where the internal boundary meets the
    // external boundary, not a manifold mesh.  We still have work to
    // do to properly support non-manifold meshes, and
    // find_neighbors() when preparing such a mesh is likely to fail,
    // so we'll just do the best test we can here, requesting a set of
    // three different boundaries.  We drop the left and bottom
    // boundaries because they have any T intersections with the
    // interior boundaries, and we drop one of the interior boundaries
    // because with both together we'd have 4 edges meeting at every
    // interior boundary vertex.
    const std::set<boundary_id_type> multi_boundaries {1, 2, 5};
    _mesh->get_boundary_info().sync(multi_boundaries,
                                    *_multi_boundary_mesh);

    // Check that the interior_parent side indices that we set on the
    // BoundaryMesh as extra integers agree with the side returned
    // by which_side_am_i().
    unsigned int parent_side_index_tag =
      _multi_boundary_mesh->get_elem_integer_index("parent_side_index");

    for (const auto & belem : _multi_boundary_mesh->element_ptr_range())
      {
        dof_id_type parent_side_index =
          belem->get_extra_integer(parent_side_index_tag);

        if (belem->interior_parent() != remote_elem)
          CPPUNIT_ASSERT_EQUAL
            (static_cast<dof_id_type>(belem->interior_parent()->which_side_am_i(belem)),
             parent_side_index);
      }

    std::set<boundary_id_type> left_id, right_id;
    left_id.insert(3);
    right_id.insert(1);

    // Add the right side of the square to the square; this should
    // make it a mixed dimension mesh. We skip storing the parent
    // side ids (which is the default) since they are not needed
    // in this particular test.
    _mesh->get_boundary_info().add_elements
      (right_id, *_mesh, /*store_parent_side_ids=*/false);
    _mesh->prepare_for_use();

    // Add the left side of the square to its own boundary mesh.
    _mesh->get_boundary_info().sync(left_id, *_left_boundary_mesh);


    // Create a BoundaryMesh from the internal sidesets relative to subdomain 1.
    {
      std::set<boundary_id_type> requested_boundary_ids {bid1, bid2};
      std::set<subdomain_id_type> subdomains_relative_to;
      subdomains_relative_to.insert(1);
      _mesh->get_boundary_info().sync(requested_boundary_ids,
                                      *_internal_boundary_mesh,
                                      subdomains_relative_to);
    }
  }

public:
  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();
    this->sync_and_test_meshes();
#endif
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    // There'd better be 3*5 + 5 elements in the interior plus right
    // boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(20),
                         _mesh->n_elem());

    // There'd better be 7*11 nodes in the interior
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(77),
                         _mesh->n_nodes());

    // There'd better be 2*(3+5) elements on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(16),
                         _exterior_boundary_mesh->n_elem());

    // There'd better be 2*2*(3+5) nodes on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(32),
                         _exterior_boundary_mesh->n_nodes());

    // There'd better be 5 elements on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(5),
                         _left_boundary_mesh->n_elem());

    // There'd better be 2*5+1 nodes on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(11),
                         _left_boundary_mesh->n_nodes());

    // There are four elements in the one-sided internal sideset mesh.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(4),
                         _internal_boundary_mesh->n_elem());

    // There are 2*n_elem + 1 nodes in the one-sided internal sideset mesh.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(9),
                         _internal_boundary_mesh->n_nodes());

    // There'd better be 3+5 + 4 elements on the multi-boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12),
                         _multi_boundary_mesh->n_elem());

    // There'd better be (3+5)*2+1 + 4*2+1 nodes on the multi-boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(26),
                         _multi_boundary_mesh->n_nodes());

    this->sanityCheck();
  }

  void testBoundarySerialization()
  {
    LOG_UNIT_TEST;

    MeshSerializer internal_serializer(*_internal_boundary_mesh);

    // There are four elements in the internal sideset mesh.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(4),
                         _internal_boundary_mesh->n_elem());

    // There are 2*n_elem + 1 nodes in the internal sideset mesh.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(9),
                         _internal_boundary_mesh->n_nodes());

    this->sanityCheck();
  }

  void sanityCheck()
  {
    // Sanity check all the elements
    for (const auto & elem : _mesh->active_element_ptr_range())
      {
        const Elem * pip = elem->dim() < 2 ? elem->interior_parent() : nullptr;

        // On a DistributedMesh we might not be able to see the
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
            LIBMESH_ASSERT_FP_EQUAL(0.8, elem->vertex_average()(0),
                                    TOLERANCE*TOLERANCE);
          }
        else
          {
            CPPUNIT_ASSERT_EQUAL(elem->type(), QUAD9);
          }
      }

    for (const auto & elem : _left_boundary_mesh->active_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);

        const Elem * pip = elem->interior_parent();

        // On a DistributedMesh we might not be able to see the
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
        LIBMESH_ASSERT_FP_EQUAL(0.2, elem->vertex_average()(0),
                                TOLERANCE*TOLERANCE);
      }

    for (const auto & elem : _left_boundary_mesh->active_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);

        const Elem * pip = elem->interior_parent();

        // On a DistributedMesh we might not be able to see the
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

    // Sanity check for the internal sideset mesh.
    for (const auto & elem : _internal_boundary_mesh->active_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->type(), EDGE3);

        // All of the elements in the internal sideset mesh should
        // have the same subdomain id as the parent Elems (i.e. 1)
        // they came from.
        CPPUNIT_ASSERT_EQUAL(static_cast<subdomain_id_type>(1),
                             elem->subdomain_id());
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
  LIBMESH_CPPUNIT_TEST_SUITE( BoundaryRefinedMeshTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testBoundarySerialization );
#endif

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();
    this->sync_and_test_meshes();

    // Need to refine interior mesh before separate boundary meshes,
    // if we want to get interior_parent links right.
    MeshRefinement(*_mesh).uniformly_refine(1);

    MeshRefinement(*_multi_boundary_mesh).uniformly_refine(1);
    MeshRefinement(*_exterior_boundary_mesh).uniformly_refine(1);
    MeshRefinement(*_left_boundary_mesh).uniformly_refine(1);
    MeshRefinement(*_internal_boundary_mesh).uniformly_refine(1);
#endif
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    // There'd better be 3*5*4 + 5*2 active elements in the interior
    // plus right boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(70),
                         _mesh->n_active_elem());

    // Plus the original 20 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(90),
                         _mesh->n_elem());

    // There'd better be 13*21 nodes in the interior
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(273),
                         _mesh->n_nodes());

    // There'd better be 2*2*(3+5) active elements on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(32),
                         _exterior_boundary_mesh->n_active_elem());

    // Plus the original 16 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(48),
                         _exterior_boundary_mesh->n_elem());

    // There'd better be 2*2*2*(3+5) nodes on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(64),
                         _exterior_boundary_mesh->n_nodes());

    // There'd better be 2*5 active elements on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(10),
                         _left_boundary_mesh->n_active_elem());

    // Plus the original 5 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(15),
                         _left_boundary_mesh->n_elem());

    // There'd better be 2*2*5+1 nodes on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(21),
                         _left_boundary_mesh->n_nodes());

    // There'd better be 2*4 active elements on the one-sided internal
    // sideset mesh
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8),
                         _internal_boundary_mesh->n_active_elem());

    // Plus the original 4 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12),
                         _internal_boundary_mesh->n_elem());

    // There'd better be 2*2*4+1 nodes on the internal boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17),
                         _internal_boundary_mesh->n_nodes());

    // There'd better be 2*(3+5) + 2*4 active elements on the
    // multi-boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(24),
                         _multi_boundary_mesh->n_active_elem());

    // Plus the original 12 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(36),
                         _multi_boundary_mesh->n_elem());

    // There'd better be 2*2*(3+5)+1 + 2*2*4+1 nodes on the
    // multi-boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(50),
                         _multi_boundary_mesh->n_nodes());


    this->sanityCheck();
  }


  void testBoundarySerialization()
  {
    LOG_UNIT_TEST;

    MeshSerializer internal_serializer(*_internal_boundary_mesh);

    // There'd better be 2*4 active elements on the internal
    // sideset mesh
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8),
                         _internal_boundary_mesh->n_active_elem());

    // Plus the original 4 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12),
                         _internal_boundary_mesh->n_elem());

    // There'd better be 2*2*4+1 nodes on the internal boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17),
                         _internal_boundary_mesh->n_nodes());

    this->sanityCheck();
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryRefinedMeshTest );


class BoundaryOfRefinedMeshTest : public BoundaryMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement *before* synchronizing to the boundary mesh
   * rather than after.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( BoundaryOfRefinedMeshTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testBoundarySerialization );
#endif

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();

    // Refine interior mesh before creating boundary meshes
    MeshRefinement(*_mesh).uniformly_refine(1);

    this->sync_and_test_meshes();
#endif
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    // There'd better be 3*5*4 + 5*2 active elements in the interior
    // plus right boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(70),
                         _mesh->n_active_elem());

    // Plus the original 20 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(90),
                         _mesh->n_elem());

    // There'd better be 13*21 nodes in the interior
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(273),
                         _mesh->n_nodes());

    // There'd better be 2*2*(3+5) active elements on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(32),
                         _exterior_boundary_mesh->n_active_elem());

    // Plus the original 16 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(48),
                         _exterior_boundary_mesh->n_elem());

    // There'd better be 2*2*2*(3+5) nodes on the exterior boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(64),
                         _exterior_boundary_mesh->n_nodes());

    // There'd better be 2*5 active elements on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(10),
                         _left_boundary_mesh->n_active_elem());

    // Plus the original 5 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(15),
                         _left_boundary_mesh->n_elem());

    // There'd better be 2*2*5+1 nodes on the left boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(21),
                         _left_boundary_mesh->n_nodes());

    // There'd better be 2*4 active elements on the internal
    // sideset mesh
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8),
                         _internal_boundary_mesh->n_active_elem());

    // Plus the original 4 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12),
                         _internal_boundary_mesh->n_elem());

    // There'd better be 2*2*4+1 nodes on the one-sided internal boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17),
                         _internal_boundary_mesh->n_nodes());

    // There'd better be 2*(3+5) + 2*4 active elements on the
    // multi-boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(24),
                         _multi_boundary_mesh->n_active_elem());

    // Plus the original 12 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(36),
                         _multi_boundary_mesh->n_elem());

    // There'd better be 2*2*2*(3+5) + 2*2*4-1 nodes on the total boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(50),
                         _multi_boundary_mesh->n_nodes());

    this->sanityCheck();
  }

  void testBoundarySerialization()
  {
    LOG_UNIT_TEST;

    MeshSerializer internal_serializer(*_internal_boundary_mesh);

    // There'd better be 2*4 active elements on the internal
    // sideset mesh
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8),
                         _internal_boundary_mesh->n_active_elem());

    // Plus the original 4 now-inactive elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12),
                         _internal_boundary_mesh->n_elem());

    // There'd better be 2*2*4+1 nodes on the internal boundary
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17),
                         _internal_boundary_mesh->n_nodes());

    this->sanityCheck();
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryOfRefinedMeshTest );

#endif
