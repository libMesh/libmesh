#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/parallel_implementation.h> // max()
#include <libmesh/simplex_refiner.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>


using namespace libMesh;


class SimplexRefinementTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * interfaces to tetrahedralization libraries
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( SimplexRefinementTest );

  CPPUNIT_TEST( testTriRefinement );
  CPPUNIT_TEST( test3DTriRefinement );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testRefinement(UnstructuredMesh & mesh,
                      Real desired_measure,
                      Real expected_total_measure)
  {
    SimplexRefiner simplex_refiner(mesh);
    simplex_refiner.desired_volume() = desired_measure;

    simplex_refiner.refine_elements();

    Real total_measure = 0;
    for (const Elem * elem : mesh.local_element_ptr_range())
      {
        Real measure = elem->volume();
        CPPUNIT_ASSERT_LESSEQUAL(desired_measure, measure);
        total_measure += measure;
      }

    TestCommWorld->sum(total_measure);
    LIBMESH_ASSERT_FP_EQUAL(total_measure, expected_total_measure,
                            TOLERANCE*TOLERANCE);
  }


  void testTriRefinement()
  {
    LOG_UNIT_TEST;

    Mesh trimesh(*TestCommWorld);

    MeshTools::Generation::build_square (trimesh, 2, 2,
                                         0.0, 2.0, 0.0, 2.0, TRI3);

    testRefinement(trimesh, 0.05, 4);
  }


  void test3DTriRefinement()
  {
    // We may have a bug when trying to edge refine on more processors
    // than elements?  I can't reproduce it and it seems intermittent
    // in CI.
    if (TestCommWorld->size() > 8)
      return;

    LOG_UNIT_TEST;

    Mesh trimesh(*TestCommWorld);

    MeshTools::Generation::surface_octahedron
      (trimesh, -1, 1, -1, 1, -1, 1);

    testRefinement(trimesh, 0.05, 4*std::sqrt(Real(3)));
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( SimplexRefinementTest );
