#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/face_quad4_shell.h>
#include <libmesh/equation_systems.h>
#include <libmesh/zero_function.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/parallel.h>
#include <libmesh/mesh_refinement.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <regex>

using namespace libMesh;

class BoundaryInfoTest : public CppUnit::TestCase {
  /**
   * This test ensures various aspects of the BoundaryInfo class work as expected.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( BoundaryInfoTest );

  CPPUNIT_TEST( testNameCopying );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testRenumber );
# ifdef LIBMESH_ENABLE_AMR
#  ifdef LIBMESH_ENABLE_EXCEPTIONS
  CPPUNIT_TEST( testBoundaryOnChildrenErrors );
#  endif
  CPPUNIT_TEST( testBoundaryOnChildrenElementsRefineCoarsen );
  CPPUNIT_TEST( testBoundaryOnChildrenBoundaryIDs );
  CPPUNIT_TEST( testBoundaryOnChildrenBoundarySides );
# endif
# ifdef LIBMESH_ENABLE_DIRICHLET
  CPPUNIT_TEST( testShellFaceConstraints );
# endif
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testEdgeBoundaryConditions );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    MeshTools::Generation::build_square(mesh,
                                        2, 2,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh.get_boundary_info();

    // Side lists should be cleared and refilled by each call
#ifdef LIBMESH_ENABLE_DEPRECATED
    std::vector<dof_id_type> element_id_list;
    std::vector<unsigned short int> side_list;
    std::vector<boundary_id_type> bc_id_list;
#endif

    // build_square adds boundary_ids 0,1,2,3 for the bottom, right,
    // top, and left sides, respectively.

    // On a ReplicatedMesh, we should see all 4 ids on each processor
    if (mesh.is_serial())
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(4), bi.n_boundary_ids());

    // On any mesh, we should see each id on *some* processor
    {
      const std::set<boundary_id_type> & bc_ids = bi.get_boundary_ids();
      for (boundary_id_type i = 0 ; i != 4; ++i)
        {
          bool has_bcid = bc_ids.count(i);
          mesh.comm().max(has_bcid);
          CPPUNIT_ASSERT(has_bcid);
        }
    }

    // Build the side list
#ifdef LIBMESH_ENABLE_DEPRECATED
    bi.build_side_list (element_id_list, side_list, bc_id_list);
#endif

    // Test that the new vector-of-tuples API works equivalently.
    auto bc_triples = bi.build_side_list();

    // Check that there are exactly 8 sides in the BoundaryInfo for a
    // replicated mesh
    if (mesh.is_serial())
      {
#ifdef LIBMESH_ENABLE_DEPRECATED
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(8), element_id_list.size());
#endif
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(8), bc_triples.size());
      }

    // Let's test that they are preserved (in a relative sense) when
    // we clone a mesh.
    std::unique_ptr<MeshBase> mesh_clone = mesh.clone();
    CPPUNIT_ASSERT(mesh_clone->get_boundary_info() ==
                   mesh.get_boundary_info());

    // Let's test that we can remove them successfully.
    bi.remove_id(0);

    CPPUNIT_ASSERT(mesh_clone->get_boundary_info() !=
                   mesh.get_boundary_info());

    // Check that there are now only 3 boundary ids total on the Mesh.
    if (mesh.is_serial())
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3), bi.n_boundary_ids());

    {
      const std::set<boundary_id_type> & bc_ids = bi.get_boundary_ids();
      CPPUNIT_ASSERT(!bc_ids.count(0));
      for (boundary_id_type i = 1 ; i != 4; ++i)
        {
          bool has_bcid = bc_ids.count(i);
          mesh.comm().max(has_bcid);
          CPPUNIT_ASSERT(has_bcid);
        }
    }

    // Build the side list again
#ifdef LIBMESH_ENABLE_DEPRECATED
    bi.build_side_list (element_id_list, side_list, bc_id_list);
#endif
    bc_triples = bi.build_side_list();

    // Check that there are now exactly 6 sides left in the
    // BoundaryInfo on a replicated mesh
    if (mesh.is_serial())
      {
#ifdef LIBMESH_ENABLE_DEPRECATED
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(6), element_id_list.size());
#endif
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(6), bc_triples.size());
      }

    // Check that the removed ID is really removed
#ifdef LIBMESH_ENABLE_DEPRECATED
    CPPUNIT_ASSERT(std::find(bc_id_list.begin(), bc_id_list.end(), 0) == bc_id_list.end());
#endif
    typedef std::tuple<dof_id_type, unsigned short int, boundary_id_type> Tuple;
    CPPUNIT_ASSERT(std::find_if(bc_triples.begin(), bc_triples.end(),
                                [](const Tuple & t)->bool { return std::get<2>(t) == 0; }) == bc_triples.end());

    // Remove the same id again, make sure nothing changes.
    bi.remove_id(0);
    if (mesh.is_serial())
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3), bi.n_boundary_ids());

    // Remove the remaining IDs, verify that we have no sides left and
    // that we can safely reuse the same vectors in the
    // build_side_list() call.
    bi.remove_id(1);
    bi.remove_id(2);
    bi.remove_id(3);
#ifdef LIBMESH_ENABLE_DEPRECATED
    bi.build_side_list (element_id_list, side_list, bc_id_list);
#endif
    bc_triples = bi.build_side_list();

    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bi.n_boundary_ids());
#ifdef LIBMESH_ENABLE_DEPRECATED
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), element_id_list.size());
#endif
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bc_triples.size());
  }


  void testRenumber()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    MeshTools::Generation::build_square(mesh,
                                        2, 2,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh.get_boundary_info();

    // Side lists should be cleared and refilled by each call
#ifdef LIBMESH_ENABLE_DEPRECATED
    std::vector<dof_id_type> element_id_list;
    std::vector<unsigned short int> side_list;
    std::vector<boundary_id_type> bc_id_list;
#endif

    // build_square adds boundary_ids 0,1,2,3 for the bottom, right,
    // top, and left sides, respectively.  Let's remap those, not 1-1.
    bi.renumber_id(0, 4);
    bi.renumber_id(1, 5);
    bi.renumber_id(2, 6);
    bi.renumber_id(3, 6);

    const std::map<boundary_id_type, std::string> expected_names =
      {{4,"bottom"}, {5,"right"}, {6,"left"}};

    // On a ReplicatedMesh, we should see ids 4,5,6 on each processor
    if (mesh.is_serial())
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3), bi.n_boundary_ids());

    // On any mesh, we should see each new id on *some* processor, and
    // shouldn't see old ids on *any* processor
    {
      const std::set<boundary_id_type> & bc_ids = bi.get_boundary_ids();
      for (boundary_id_type i = 0 ; i != 4; ++i)
        {
          bool has_bcid = bc_ids.count(i);
          mesh.comm().max(has_bcid);
          CPPUNIT_ASSERT(!has_bcid);
        }
      for (boundary_id_type i = 4 ; i != 7; ++i)
        {
          bool has_bcid = bc_ids.count(i);

          bool bad_name = false;
          if (has_bcid)
          {
            const std::string & current_name = bi.sideset_name(i);

            bad_name = (current_name != libmesh_map_find(expected_names, i));
          }

          // At least one proc should have each of these BCs
          mesh.comm().max(has_bcid);
          CPPUNIT_ASSERT(has_bcid);

          // No proc should have the wrong name for a BC it has
          mesh.comm().max(bad_name);
          CPPUNIT_ASSERT(!bad_name);
        }
    }

    // Check that there are still exactly 8 sides in the BoundaryInfo
    // for a replicated mesh
    auto bc_triples = bi.build_side_list();

    if (mesh.is_serial())
      {
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(8), bc_triples.size());
      }

    // Remove the new IDs, verify that we have no sides left
    bi.remove_id(4);
    bi.remove_id(5);
    bi.remove_id(6);
    bc_triples = bi.build_side_list();

    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bi.n_boundary_ids());
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bc_triples.size());
  }


  void testEdgeBoundaryConditions()
  {
    LOG_UNIT_TEST;

    const unsigned int n_elem = 5;
    const std::string mesh_filename = "cube_mesh.xda";

    {
      Mesh mesh(*TestCommWorld);
      MeshTools::Generation::build_cube(mesh,
                                        n_elem, n_elem, n_elem,
                                        0., 1.,
                                        0., 1.,
                                        0., 1.,
                                        HEX8);

      BoundaryInfo & bi = mesh.get_boundary_info();

      // build_cube does not add any edge boundary IDs
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bi.n_edge_conds());

      // Let's now add some edge boundary IDs.
      // We loop over all elements (not just local elements) so that
      // all processors know about the boundary IDs
      const boundary_id_type BOUNDARY_ID_MAX_X = 2;
      const boundary_id_type BOUNDARY_ID_MIN_Y = 1;
      const boundary_id_type EDGE_BOUNDARY_ID = 20;

      for (const auto & elem : mesh.element_ptr_range())
        {
          unsigned short side_max_x = 0, side_min_y = 0;
          bool found_side_max_x = false, found_side_min_y = false;

          for (unsigned short side=0; side<elem->n_sides(); side++)
            {
              if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
                {
                  side_max_x = side;
                  found_side_max_x = true;
                }

              if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MIN_Y))
                {
                  side_min_y = side;
                  found_side_min_y = true;
                }
            }

          // If elem has sides on boundaries
          // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
          // then let's set an edge boundary condition
          if (found_side_max_x && found_side_min_y)
            for (unsigned short e=0; e<elem->n_edges(); e++)
              if (elem->is_edge_on_side(e, side_max_x) &&
                  elem->is_edge_on_side(e, side_min_y))
                bi.add_edge(elem, e, EDGE_BOUNDARY_ID);
        }

      // Check that we have the expected number of edge boundary IDs after
      // updating bi
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(n_elem), bi.n_edge_conds());

      // Let's test that edge BCIDs are preserved (in a relative
      // sense) when we clone a mesh.
      std::unique_ptr<MeshBase> mesh_clone = mesh.clone();
      CPPUNIT_ASSERT(mesh_clone->get_boundary_info() ==
                     mesh.get_boundary_info());

      mesh.write(mesh_filename);
    }

    // Make sure all processors are done writing before we try to
    // start reading
    TestCommWorld->barrier();

    Mesh mesh(*TestCommWorld);
    mesh.read(mesh_filename);

    // Check that writing and reading preserves the edge boundary IDs
    BoundaryInfo & bi = mesh.get_boundary_info();
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(n_elem), bi.n_edge_conds());
  }

  void testNameCopying()
  {
    LOG_UNIT_TEST;

      Mesh mesh(*TestCommWorld);
      MeshTools::Generation::build_line(mesh,
                                        8,
                                        0., 1.,
                                        EDGE2);

      Mesh mesh2(mesh);

      BoundaryInfo & bi = mesh.get_boundary_info();
      bi.sideset_name(0) = "zero";
      bi.sideset_name(1) = "one";
      bi.sideset_name(2) = "two";
      bi.sideset_name(3) = "three";
      bi.nodeset_name(0) = "ZERO";
      bi.nodeset_name(1) = "ONE";

      BoundaryInfo bi2 {bi};
      CPPUNIT_ASSERT_EQUAL(bi2.get_sideset_name(0), std::string("zero"));
      CPPUNIT_ASSERT_EQUAL(bi2.get_sideset_name(1), std::string("one"));
      CPPUNIT_ASSERT_EQUAL(bi2.get_sideset_name(2), std::string("two"));
      CPPUNIT_ASSERT_EQUAL(bi2.get_sideset_name(3), std::string("three"));
      CPPUNIT_ASSERT_EQUAL(bi2.get_nodeset_name(0), std::string("ZERO"));
      CPPUNIT_ASSERT_EQUAL(bi2.get_nodeset_name(1), std::string("ONE"));

      BoundaryInfo & bi3 = mesh2.get_boundary_info();
      bi3 = bi;
      CPPUNIT_ASSERT_EQUAL(bi3.get_sideset_name(0), std::string("zero"));
      CPPUNIT_ASSERT_EQUAL(bi3.get_sideset_name(1), std::string("one"));
      CPPUNIT_ASSERT_EQUAL(bi3.get_sideset_name(2), std::string("two"));
      CPPUNIT_ASSERT_EQUAL(bi3.get_sideset_name(3), std::string("three"));
      CPPUNIT_ASSERT_EQUAL(bi3.get_nodeset_name(0), std::string("ZERO"));
      CPPUNIT_ASSERT_EQUAL(bi3.get_nodeset_name(1), std::string("ONE"));
  }

#ifdef LIBMESH_ENABLE_DIRICHLET
  void testShellFaceConstraints()
  {
    LOG_UNIT_TEST;

    // Make a simple two element mesh that we can use to test constraints
    Mesh mesh(*TestCommWorld);

    //  (0,1)           (1,1)
    //  x---------------x
    //  |               |
    //  |               |
    //  |               |
    //  |               |
    //  |               |
    //  x---------------x
    //  (0,0)           (1,0)
    //  |               |
    //  |               |
    //  |               |
    //  |               |
    //  x---------------x
    //  (0,-1)          (1,-1)

    mesh.add_point( Point(0.0,-1.0), 4 );
    mesh.add_point( Point(1.0,-1.0), 5 );
    mesh.add_point( Point(1.0, 0.0), 1 );
    mesh.add_point( Point(1.0, 1.0), 2 );
    mesh.add_point( Point(0.0, 1.0), 3 );
    mesh.add_point( Point(0.0, 0.0), 0 );

    Elem * elem_top = mesh.add_elem(Elem::build(QUADSHELL4));
    elem_top->set_node(0) = mesh.node_ptr(0);
    elem_top->set_node(1) = mesh.node_ptr(1);
    elem_top->set_node(2) = mesh.node_ptr(2);
    elem_top->set_node(3) = mesh.node_ptr(3);

    Elem * elem_bottom = mesh.add_elem(Elem::build(QUADSHELL4));
    elem_bottom->set_node(0) = mesh.node_ptr(4);
    elem_bottom->set_node(1) = mesh.node_ptr(5);
    elem_bottom->set_node(2) = mesh.node_ptr(1);
    elem_bottom->set_node(3) = mesh.node_ptr(0);

    BoundaryInfo & bi = mesh.get_boundary_info();
    bi.add_shellface(elem_top, 0, 10);
    bi.add_shellface(elem_bottom, 1, 20);

    mesh.allow_renumbering(true);
    mesh.prepare_for_use();

    // Let's test that shellface BCIDs are preserved (in a relative
    // sense) when we clone a mesh.
    std::unique_ptr<MeshBase> mesh_clone = mesh.clone();
    CPPUNIT_ASSERT(mesh_clone->get_boundary_info() ==
                   mesh.get_boundary_info());

    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(2), bi.n_shellface_conds());

    EquationSystems es(mesh);
    System & system = es.add_system<System> ("SimpleSystem");
    system.add_variable("u", FIRST);

    // Add a Dirichlet constraint to check that we impose constraints
    // correctly on shell faces.
    std::vector<unsigned int> variables;
    variables.push_back(0);
    std::set<boundary_id_type> shellface_ids;
    shellface_ids.insert(20);
    ZeroFunction<> zf;
    DirichletBoundary dirichlet_bc(shellface_ids,
                                   variables,
                                   &zf);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    es.init();

    // Find elem_bottom again if we have it (it may have been deleted
    // in a DistributedMesh or renumbered in theory)
    elem_bottom = nullptr;
    for (unsigned int e = 0; e != mesh.max_elem_id(); ++e)
      {
        Elem *elem = mesh.query_elem_ptr(e);
        if (elem && elem->point(3) == Point(0,0))
          elem_bottom = elem;
      }

    // We expect to have a dof constraint on all four dofs of
    // elem_bottom
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(4), static_cast<std::size_t>(system.n_constrained_dofs()));

    // But we may only know the details of that
    // constraint on the processor which owns elem_bottom.
    if (elem_bottom &&
        elem_bottom->processor_id() == mesh.processor_id())
      {
        std::vector<dof_id_type> dof_indices;
        system.get_dof_map().dof_indices(elem_bottom, dof_indices);
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(4), dof_indices.size());

        for(unsigned int i=0; i<dof_indices.size(); i++)
          {
            dof_id_type dof_id = dof_indices[i];
            CPPUNIT_ASSERT( system.get_dof_map().is_constrained_dof(dof_id) );
          }
      }
  }
#endif // LIBMESH_ENABLE_DIRICHLET

#if LIBMESH_ENABLE_AMR
#  if LIBMESH_ENABLE_EXCEPTIONS
  void testBoundaryOnChildrenErrors()
  {
    LOG_UNIT_TEST;

    // We create one cell only. The default boundaries of the cell are below.
    //   ___2___
    // 3 |     | 1
    //   |_____|
    //      0

    auto mesh = std::make_unique<Mesh>(*TestCommWorld);
    MeshTools::Generation::build_square(*mesh,
                                        1, 1,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh->get_boundary_info();

    // We only have one element, but for easy access we use the iterator
    for (auto & elem : mesh->active_element_ptr_range())
      elem->set_refinement_flag(Elem::REFINE);
    mesh->prepare_for_use();

    MeshRefinement(*mesh).refine_elements();
    mesh->prepare_for_use();

    // Now we try to add boundary id 3 to a child on side 3. This should
    // result in a "not implemented" error message
    bool threw_desired_exception = false;
    try {
      for (auto & elem : mesh->active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        if (c(0) < 0.5 && c(1) > 0.5)
          bi.add_side(elem, 3, 3);
      }
    }
    catch (libMesh::NotImplemented & e) {
      std::regex msg_regex("Trying to add boundary ID 3 which already exists on the ancestors");
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_desired_exception = true;
    }
    // If we have more than 4 processors, or a poor partitioner, we
    // might not get an exception on every processor
    mesh->comm().max(threw_desired_exception);

    CPPUNIT_ASSERT(threw_desired_exception);

    threw_desired_exception = false;
    try {
      for (auto & elem : mesh->active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        if (c(0) < 0.5 && c(1) > 0.5)
          bi.add_side(elem, 3, {3,4});
      }
    }
    catch (libMesh::NotImplemented & e) {
      std::regex msg_regex("Trying to add boundary ID 3 which already exists on the ancestors");
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_desired_exception = true;
    }

    // If we have more than 4 processors, or a poor partitioner, we
    // might not get an exception on every processor
    mesh->comm().max(threw_desired_exception);

    CPPUNIT_ASSERT(threw_desired_exception);

    // We tested the side addition errors, now we move to the removal parts.
    // We will attempt the removal of boundary 3 through the child
    threw_desired_exception = false;
    bi.allow_children_on_boundary_side(true);
    try {
      for (auto & elem : mesh->active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        if (c(0) < 0.5 && c(1) > 0.5)
          bi.remove_side(elem, 3, 3);
      }
    }
    catch (libMesh::NotImplemented & e) {
      std::regex msg_regex("We cannot delete boundary ID 3 using a child because it is inherited from an ancestor");
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_desired_exception = true;
    }

    // If we have more than 4 processors, or a poor partitioner, we
    // might not get an exception on every processor
    mesh->comm().max(threw_desired_exception);

    CPPUNIT_ASSERT(threw_desired_exception);
  }
#  endif // LIBMESH_ENABLE_EXCEPTIONS

  void testBoundaryOnChildrenElementsRefineCoarsen()
  {
    LOG_UNIT_TEST;

    // Set subdomain ids for specific elements, we will refine/coarsen
    // the cell on subdomain 1
    // _____________
    // |  1  |  2  |
    // |_____|_____|

    auto mesh = std::make_unique<Mesh>(*TestCommWorld);
    MeshTools::Generation::build_square(*mesh,
                                        2, 1,
                                        0., 2.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh->get_boundary_info();

    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1)
      {
        elem->subdomain_id() = 1;
        elem->set_refinement_flag(Elem::REFINE);
      }
      else
        elem->subdomain_id() = 2;
    }
    mesh->prepare_for_use();

    // Refine the elements once in subdomain 1, and
    // add the right side subdomain 1 as boundary 5
    MeshRefinement(*mesh).refine_elements();
    mesh->prepare_for_use();

    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1 && c(0) > 0.5)
        bi.add_side(elem, 1, 5);
    }
    mesh->prepare_for_use();

    // Check the middle boundary, we expect to have two sides in boundary 5
    unsigned int count = 0;
    for (auto & elem : mesh->active_element_ptr_range())
      if (bi.has_boundary_id(elem, 1, 5))
        count++;

    if (mesh->n_active_local_elem())
    {
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), count);
      CPPUNIT_ASSERT(bi.is_children_on_boundary_side());
    }

    // First, we will coarsen the the elements on subdomain 1. This
    // is to check if the boundary information propagates upward upon
    // coarsening.
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1)
        elem->set_refinement_flag(Elem::COARSEN);
    }
    mesh->prepare_for_use();

    // The coarsened element should have its side on boundary 5
    // This is boundary info transferred from this child element
    MeshRefinement(*mesh).coarsen_elements();
    mesh->prepare_for_use();

    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1)
      {
        CPPUNIT_ASSERT(bi.has_boundary_id(elem, 1, 5));
        // We clean up this boundary ID for the next round of tests
        bi.remove_side(elem, 1, 5);
        // we will refine this element again
        elem->set_refinement_flag(Elem::REFINE);
      }
    }

    MeshRefinement(*mesh).refine_elements();
    mesh->prepare_for_use();

    // This time we remove boundary 5 from one of the children. We expect
    // the boundary not to propagate to the next level. Furthermore we
    // expect boundary 5 to be deleted from the parent's boundaries
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1)
        elem->set_refinement_flag(Elem::COARSEN);
      if (c(0) > 0.5 && c(0) < 1 && c(1) < 0.5)
        bi.add_side(elem, 1, 5);
    }
    mesh->prepare_for_use();

    MeshRefinement(*mesh).coarsen_elements();

    mesh->prepare_for_use();

    // The parent element should not have any side associated with boundary 5
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 1)
        CPPUNIT_ASSERT(!bi.has_boundary_id(elem, 1, 5));
    }
  }

  void testBoundaryOnChildrenBoundaryIDs()
  {
    LOG_UNIT_TEST;

    // We create one cell only. The default boundaries of the cell are below.
    // We will refine the mesh and add a new boundary id to the left side (side 3).
    // Then will query the available boundary ids on the added side. It should return
    // both the parent's and the child's boundaries.
    //   ___2___
    // 3 |     | 1
    //   |_____|
    //      0

    auto mesh = std::make_unique<Mesh>(*TestCommWorld);
    MeshTools::Generation::build_square(*mesh,
                                        1, 1,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh->get_boundary_info();

    // We only have one element, but for easy access we use the iterator
    for (auto & elem : mesh->active_element_ptr_range())
      elem->set_refinement_flag(Elem::REFINE);
    mesh->prepare_for_use();

    MeshRefinement(*mesh).refine_elements();

    // Now we add the extra boundary ID (5) to the element in the top
    // left corner
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(1) > 0.5)
        bi.add_side(elem, 3, 5);
    }
    mesh->prepare_for_use();

    // Okay, now we query the boundary ids on side 3 of the child and check if it has
    // the right elements
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(1) > 0.5)
      {
        std::vector<boundary_id_type> container;
        bi.boundary_ids(elem, 3, container);

        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(2), container.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<boundary_id_type>(5), container[0]);
        CPPUNIT_ASSERT_EQUAL(static_cast<boundary_id_type>(3), container[1]);
      }
    }
  }

  void testBoundaryOnChildrenBoundarySides()
  {
    LOG_UNIT_TEST;

    // We create one cell only. The default boundaries of the cell are below.
    // We will refine mesh and see if we can get back the correct sides
    // for a given boundary id on an internal boundary.
    //   ___2___
    // 3 |     | 1
    //   |_____|
    //      0

    auto mesh = std::make_unique<Mesh>(*TestCommWorld);
    MeshTools::Generation::build_square(*mesh,
                                        1, 1,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);

    BoundaryInfo & bi = mesh->get_boundary_info();


    // We only have one element, but for easy access we use the iterator
    for (auto & elem : mesh->active_element_ptr_range())
      elem->set_refinement_flag(Elem::REFINE);
    mesh->prepare_for_use();
    MeshRefinement(*mesh).refine_elements();

    // Now we add the extra boundary ID (5) to two sides of
    // the element in the bottom left corner. then we refine again
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(1) < 0.5)
      {
        bi.add_side(elem, 1, 5);
        bi.add_side(elem, 2, 5);
        elem->set_refinement_flag(Elem::REFINE);
      }
    }
    mesh->prepare_for_use();
    MeshRefinement(*mesh).refine_elements();

    // Okay, now we add another boundary id (6) to the cell which is in the bottom
    // right corner of the refined element
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(0) > 0.25 && c(1) < 0.25)
        bi.add_side(elem, 1, 6);
    }

    // Time to test if we can get back the boundary sides, first we
    // check if we can get back boundary from the ancestors of (5) on
    // the cell which only has boundary (6) registered. We also check
    // if we can get boundary (6) back.

    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(0) > 0.25 && c(1) < 0.25)
      {
        const auto side_5 = bi.side_with_boundary_id(elem, 5);
        const auto side_6 = bi.side_with_boundary_id(elem, 6);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), side_5);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), side_6);
      }
    }

    // Now we go and try to query the sides with boundary id (5) using
    // the element which is at the top right corner of the bottom
    // right parent.
    for (auto & elem : mesh->active_element_ptr_range())
    {
      const Point c = elem->vertex_average();
      if (c(0) < 0.5 && c(0) > 0.25 && c(1) > 0.25 && c(1) < 0.5)
      {
        const auto sides = bi.sides_with_boundary_id(elem, 5);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned long>(2), sides.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1), sides[0]);
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(2), sides[1]);
      }
    }
  }
#endif //LIBMESH_ENABLE_AMR
};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryInfoTest );
