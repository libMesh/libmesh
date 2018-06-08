// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/exodusII_io.h>

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


Number x_plus_y (const Point& p,
                 const Parameters&,
                 const std::string&,
                 const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);

  return x + y;
}


class MeshInputTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( MeshInputTest );

#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST( testExodusCopyElementSolution );
#endif

  CPPUNIT_TEST( testMeshMoveConstructor );
  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_EXODUS_API
  void testExodusCopyElementSolution ()
  {
    {
      Mesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("e", CONSTANT, MONOMIAL);

      MeshTools::Generation::build_square (mesh,
                                           3, 3,
                                           0., 1., 0., 1.);

      es.init();
      sys.project_solution(x_plus_y, NULL, es.parameters);

      ExodusII_IO exii(mesh);

      // Don't try to write element data as nodal data
      std::set<std::string> sys_list;
      exii.write_equation_systems("mesh_with_soln.e", es, &sys_list);

      // Just write it as element data
      exii.write_element_data(es);
    }

    // copy_elemental_solution currently requires ReplicatedMesh
    {
      ReplicatedMesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("teste", CONSTANT, MONOMIAL);

      ExodusII_IO exii(mesh);

      if (mesh.processor_id() == 0)
        exii.read("mesh_with_soln.e");
      MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      // Read the solution e into variable teste.
      //
      // With complex numbers, we'll only bother reading the real
      // part.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      exii.copy_elemental_solution(sys, "teste", "r_e");
#else
      exii.copy_elemental_solution(sys, "teste", "e");
#endif

      for (Real x = 1.L/6.L; x < 1; x += 1.L/3.L)
        for (Real y = 1.L/6.L; y < 1; y += 1.L/3.L)
          {
            Point p(x,y);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p)),
                                         libmesh_real(x+y),
                                         TOLERANCE*TOLERANCE);
          }
    }
  }
#endif

  void testMeshMoveConstructor ()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square (mesh,
                                         3, 3,
                                         0., 1., 0., 1.);

    // Construct mesh2, stealing the resources of the original.
    Mesh mesh2(std::move(mesh));

    // Make sure mesh2 now has the 9 elements.
    CPPUNIT_ASSERT_EQUAL(mesh2.n_elem(),
                         static_cast<unsigned int>(9));

    // Verify that the moved-from mesh's Partitioner and BoundaryInfo
    // objects were successfully stolen.  Note: moved-from unique_ptrs
    // are guaranteed to compare equal to nullptr, see e.g. Section
    // 20.8.1/4 of the standard.
    // https://stackoverflow.com/questions/24061767/is-unique-ptr-guaranteed-to-store-nullptr-after-move
    CPPUNIT_ASSERT(!mesh.partitioner());
    CPPUNIT_ASSERT(!mesh.boundary_info);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshInputTest );
