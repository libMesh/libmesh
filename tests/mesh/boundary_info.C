// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>

#include "test_comm.h"

using namespace libMesh;

class BoundaryInfoTest : public CppUnit::TestCase {
  /**
   * This test ensures various aspects of the BoundaryInfo class work as expected.
   */
public:
  CPPUNIT_TEST_SUITE( BoundaryInfoTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testEdgeBoundaryConditions );

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

    // build_square adds boundary_ids 0,1,2,3 for the bottom, right,
    // top, and left sides, respectively.  Let's test that we can
    // remove them successfully.
    BoundaryInfo & bi = mesh.get_boundary_info();
    bi.remove_id(0);

    // Check that there are now only 3 boundary ids total on the Mesh.
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3), bi.n_boundary_ids());

    // Build the side list
    std::vector<dof_id_type> element_id_list;
    std::vector<unsigned short int> side_list;
    std::vector<boundary_id_type> bc_id_list;
    bi.build_side_list (element_id_list, side_list, bc_id_list);

    // Check that there are now exactly 6 sides left in the BoundaryInfo
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(6), element_id_list.size());

    // Remove the same id again, make sure nothing changes.
    bi.remove_id(0);
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3), bi.n_boundary_ids());

    // Remove the remaining IDs, verify that we have no sides left and
    // that we can safely reuse the same vectors in the
    // build_side_list() call.
    bi.remove_id(1);
    bi.remove_id(2);
    bi.remove_id(3);
    bi.build_side_list (element_id_list, side_list, bc_id_list);

    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), bi.n_boundary_ids());
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(0), element_id_list.size());
  }

  void testEdgeBoundaryConditions()
  {
    const unsigned int n_elem = 5;
    const std::string mesh_filename = "cube_mesh.xdr";

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

      MeshBase::const_element_iterator       el     = mesh.elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.elements_end();
      for ( ; el != end_el; ++el)
        {
          const Elem * elem = *el;

          unsigned int side_max_x = 0, side_min_y = 0;
          bool found_side_max_x = false, found_side_min_y = false;

          for (unsigned int side=0; side<elem->n_sides(); side++)
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
            for (unsigned int e=0; e<elem->n_edges(); e++)
              if (elem->is_edge_on_side(e, side_max_x) &&
                  elem->is_edge_on_side(e, side_min_y))
                bi.add_edge(elem, e, EDGE_BOUNDARY_ID);
        }

      // Check that we have the expected number of edge boundary IDs after
      // updating bi
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(n_elem), bi.n_edge_conds());

      mesh.write(mesh_filename);
    }

    Mesh mesh(*TestCommWorld);
    mesh.read(mesh_filename);

    // Check that writing and reading preserves the edge boundary IDs
    BoundaryInfo & bi = mesh.get_boundary_info();
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(n_elem), bi.n_edge_conds());
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryInfoTest );
