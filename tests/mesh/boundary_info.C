// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"

using namespace libMesh;

class BoundaryInfoTest : public CppUnit::TestCase {
  /**
   * This test ensures various aspects of the BoundaryInfo class work as expected.
   */
public:
  CPPUNIT_TEST_SUITE( BoundaryInfoTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

protected:

  Mesh* _mesh;

  void build_mesh()
  {
    _mesh = new Mesh(*TestCommWorld);

    MeshTools::Generation::build_square(*_mesh,
                                        2, 2,
                                        0., 1.,
                                        0., 1.,
                                        QUAD4);
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
    // build_square adds boundary_ids 0,1,2,3 for the bottom, right,
    // top, and left sides, respectively.  Let's test that we can
    // remove them successfully.
    BoundaryInfo & bi = _mesh->get_boundary_info();
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
};

CPPUNIT_TEST_SUITE_REGISTRATION( BoundaryInfoTest );
