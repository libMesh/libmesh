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

// The systems and solvers we may use
#include "libmesh/diff_solver.h"
#include "libmesh/steady_solver.h"

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
        mesh, 2, 2, -1., 1., -1., 1., Utility::string_to_enum<ElemType>("TRI6"));

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    ExplicitSystem & system = equation_systems.add_system<ExplicitSystem>("Tester");
    system.add_variable("u", FIRST, NEDELEC_ONE);
    system.add_variable("v", FIRST, LAGRANGE);

    // Initialize the system
    equation_systems.init();

    std::vector<Real> gold_solution({0.000000000010418385555305976,
                                     0.1743561867338195,
                                     -0.041899936369794721,
                                     -0.19643153680419781,
                                     -0.0000000000067936995445650305,
                                     -0.000000000010194851982520449,
                                     -0.000000000010194851982520419,
                                     0.20284410873063882,
                                     0.19643153680181882,
                                     -0.17435618673478825,
                                     -0.13942148221042558,
                                     0.0000000000087210975169893594,
                                     0.0000000000065418642410839133,
                                     0.000000000010446687286090026,
                                     -0.021522690153508699,
                                     -0.0000000000089446310895455743,
                                     0,
                                     0,
                                     0.25000000000000006,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0});

    NumericVector<Number> & sys_solution = *(system.solution);
    for (unsigned i = 0; i < gold_solution.size(); i++)
      sys_solution.set(i, gold_solution[i]);

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

    std::vector<Real> gold_vector(
        {-0.020949, -0.020949, 0., 0.009495,  0.183852,  0.,       -0.010570, -0.010570, 0.250000,
         0.000000,  0.066228,  0., 0.002165,  0.174356,  0.125000, -0.113646, 0.071746,  0.125000,
         -0.177428, 0.019003,  0., -0.196431, 0.009296,  0.125000, -0.119165, -0.000000, 0.,
         -0.000000, -0.000000, 0., 0.198569,  0.002137,  0.,       -0.000000, 0.101422,  0.,
         0.101422,  -0.000000, 0., 0.106940,  0.095903,  0.,       0.196431,  -0.009296, 0.125000,
         0.011644,  -0.162711, 0., -0.002165, -0.174356, 0.125000, -0.075229, -0.064191, 0.,
         0.000000,  0.000000,  0., 0.000000,  -0.069710, 0.,       -0.069710, 0.000000,  0.,
         -0.010761, -0.010761, 0., 0.087454,  0.000000,  0.,       0.081935,  -0.103458, 0.125000,
         -0.000000, -0.097939, 0});

    Real tol = 1e-5;
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
