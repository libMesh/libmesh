// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/dof_map.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/mesh_refinement.h>

#include "test_comm.h"

using namespace libMesh;

class MeshRefineTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( MeshRefineTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( uniformlyRefine );
  CPPUNIT_TEST( testRefinedMesh );

  CPPUNIT_TEST_SUITE_END();

private:

  SerialMesh* _mesh;

public:
  void setUp()
  {
     /*
      (0,1)           (1,1)
        x---------------x
        |               |
        |               |
        |               |
        |               |
        |               |
        x---------------x
       (0,0)           (1,0)
     */

    _mesh = new SerialMesh(*TestCommWorld);
    MeshTools::Generation::build_square(*_mesh, 1, 1 );
  }

  void tearDown()
  {
    delete _mesh;
   }

  void testMesh()
  {
    // We should have 1 total and 1 active elements now.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)1, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)1, _mesh->n_active_elem() );
  }

  void uniformlyRefine()
  {
    /*
      (0,1)           (1,1)
        x-------x-------x
        |       |       |
        |       |       |
        x-------x-------x
        |       |       |
        |       |       |
        x-------x-------x
      (0,0)           (1,0)
     */

    MeshRefinement(*_mesh).uniformly_refine(1);
    _mesh->print_info();

    // We should have 5 total and 4 active elements now.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)5, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)4, _mesh->n_active_elem() );

  }

  void testRefinedMesh()
  {
    _mesh->print_info();

    // We should *STILL* have 13 total and 10 active elements now.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)5, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)4, _mesh->n_active_elem() );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshRefineTest );
