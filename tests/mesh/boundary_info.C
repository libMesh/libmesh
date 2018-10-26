// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/face_quad4_shell.h>
#include <libmesh/equation_systems.h>
#include <libmesh/zero_function.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>

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

class BoundaryInfoTest : public CppUnit::TestCase {
  /**
   * This test ensures various aspects of the BoundaryInfo class work as expected.
   */
public:
  CPPUNIT_TEST_SUITE( BoundaryInfoTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testShellFaceConstraints );
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

    // Let's test that we can remove them successfully.
    bi.remove_id(0);

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

  void testEdgeBoundaryConditions()
  {
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

  void testShellFaceConstraints()
  {
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

    Elem* elem_top = mesh.add_elem( new QuadShell4 );
    elem_top->set_node(0) = mesh.node_ptr(0);
    elem_top->set_node(1) = mesh.node_ptr(1);
    elem_top->set_node(2) = mesh.node_ptr(2);
    elem_top->set_node(3) = mesh.node_ptr(3);

    Elem* elem_bottom = mesh.add_elem( new QuadShell4 );
    elem_bottom->set_node(0) = mesh.node_ptr(4);
    elem_bottom->set_node(1) = mesh.node_ptr(5);
    elem_bottom->set_node(2) = mesh.node_ptr(1);
    elem_bottom->set_node(3) = mesh.node_ptr(0);

    BoundaryInfo & bi = mesh.get_boundary_info();
    bi.add_shellface(elem_top, 0, 10);
    bi.add_shellface(elem_bottom, 1, 20);

    mesh.prepare_for_use(false /*skip_renumber*/);

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

};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryInfoTest );
