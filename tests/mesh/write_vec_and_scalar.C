#include "libmesh/exodusII_io.h"

// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/numeric_vector.h"

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

// Bring in everything from the libMesh namespace
using namespace libMesh;

class WriteVecAndScalar : public CppUnit::TestCase
{
  /**
   * This test ensures you can write both vector and scalar variables
   */
public:
  CPPUNIT_TEST_SUITE(WriteVecAndScalar);

  CPPUNIT_TEST(testWrite);

  CPPUNIT_TEST_SUITE_END();

  void testWrite()
  {
    Mesh mesh(*TestCommWorld);

    MeshTools::Generation::build_square(
        mesh, 1, 1, -1., 1., -1., 1., Utility::string_to_enum<ElemType>("TRI6"));

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    ExplicitSystem & system = equation_systems.add_system<ExplicitSystem>("Tester");
    system.add_variable("u", FIRST, NEDELEC_ONE);
    system.add_variable("v", FIRST, LAGRANGE);

    // Initialize the system
    equation_systems.init();

    NumericVector<Number> & sys_solution = *(system.solution);
    for (unsigned i = 0; i < 9; i++)
      sys_solution.set(i, Real(2 * i));

#ifdef LIBMESH_HAVE_EXODUS_API

    // We write the file in the ExodusII format.
    ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);
    ExodusII_IO exio(read_mesh);
    exio.read("out.e");
    read_mesh.prepare_for_use();
    EquationSystems es2(read_mesh);
    ExplicitSystem & sys2 = es2.add_system<ExplicitSystem>("Test");

    const std::vector<std::string> & nodal_vars(exio.get_nodal_var_names());

    for (const auto & var_name : nodal_vars)
      sys2.add_variable(var_name, SECOND, LAGRANGE);

    es2.init();

    for (const auto & var_name : nodal_vars)
      exio.copy_nodal_solution(sys2, var_name, var_name, 1);

    std::vector<Real> gold_vector({-1, 3,   10,  0,  1, 12, 2,  0, 14,  0,  1.5, 11, 0.5, 1,
                                   13, 0.5, 1.5, 12, 3, 4,  16, 3, 1.5, 15, 0.5, 4,  13});

    Real tol = 1e-12;
    NumericVector<Number> & sys2_soln(*sys2.solution);
    for (unsigned i = 0; i < gold_vector.size(); i++)
    {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(3 * i)), gold_vector[i], tol);
#else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(i)), gold_vector[i], tol);
#endif // #ifdef LIBMESH_USE_COMPLEX_NUMBERS
    }

#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteVecAndScalar);
