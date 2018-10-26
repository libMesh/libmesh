// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include "test_comm.h"

// C++ includes
#include <unordered_set>
#include <unordered_map>

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

class DistortTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to make sure that boundary nodes are not
   * restricted during distortion.
   */
public:
  CPPUNIT_TEST_SUITE( DistortTest );

  // 2D tests
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testDistortQuad );
#endif

  // 3D tests
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testDistortHex );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  // Helper function called by the test implementations, saves a few lines of code.
  void test_helper_2D(ElemType elem_type)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/2, /*ny=*/2,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        elem_type);

    perturb_and_check(mesh);
  }

  // Helper function called by the test implementations in 3D, saves a few lines of code.
  void test_helper_3D(ElemType elem_type)
  {
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/3);

    MeshTools::Generation::build_cube(mesh,
                                      /*nx=*/2, /*ny=*/2, /*nz=*/2,
                                      /*xmin=*/0., /*xmax=*/1.,
                                      /*ymin=*/0., /*ymax=*/1.,
                                      /*zmin=*/0., /*zmax=*/1.,
                                      elem_type);
    perturb_and_check(mesh);
  }

  // Code to save node positions, perturb the mesh, and ensure the correct nodes moved.
  void perturb_and_check(ReplicatedMesh & mesh)
  {
    // Record node positions ahead of time to make sure the correct
    // ones are moved by distort.
    std::unordered_map<dof_id_type, Point> pts_before;
    for (const auto & node : mesh.node_ptr_range())
      pts_before[node->id()] = *node;

    std::unordered_set<dof_id_type> boundary_node_ids =
      MeshTools::find_boundary_nodes (mesh);

    MeshTools::Modification::distort(mesh,
                                     /*factor=*/0.1,
                                     /*perturb_boundary=*/false);

    // Make sure the boundary nodes are not perturbed, and the
    // other nodes are.
    for (const auto & node : mesh.node_ptr_range())
      {
        bool equal = node->absolute_fuzzy_equals(pts_before[node->id()]);
        CPPUNIT_ASSERT(boundary_node_ids.count(node->id()) ? equal : !equal);
      }
  }

public:
  void setUp() {}
  void tearDown() {}
  void testDistortQuad() { test_helper_2D(QUAD4); }
  void testDistortHex() { test_helper_3D(HEX8); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( DistortTest );
