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
#include "libmesh/node.h"
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

    // We set our initial conditions based on build_square node ids
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_square(mesh,
                                        1, 1,
                                        -1., 1.,
                                        -1., 1.,
                                        Utility::string_to_enum<ElemType>("TRI6"));

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    ExplicitSystem & system = equation_systems.add_system<ExplicitSystem>("Tester");
    system.add_variable("u", FIRST, NEDELEC_ONE);
    system.add_variable("v", FIRST, LAGRANGE);

    // Initialize the system
    equation_systems.init();

    // Mimic the initial conditions of u(i)=2*i in the serial case,
    // but indexed by node id rather than dof id.
    std::vector<Real> initial_vector({10, 0, 12, 8, 4, 2, 16, 6, 14});

    NumericVector<Number> & sys_solution = *(system.solution);
    for (const auto & node : mesh.node_ptr_range())
      {
        dof_id_type dof_id;
        if (node->n_comp(0,0))
          {
            dof_id = node->dof_number(0,0,0);
          }
        else
          {
            CPPUNIT_ASSERT(node->n_comp(0,1));
            dof_id = node->dof_number(0,1,0);
          }
        if (dof_id >= sys_solution.first_local_index() &&
            dof_id < sys_solution.last_local_index())
          sys_solution.set(dof_id, initial_vector[node->id()]);
        if (mesh.n_processors() == 1)
          CPPUNIT_ASSERT_EQUAL(initial_vector[node->id()], dof_id*Real(2));
      }
    sys_solution.close();

#ifdef LIBMESH_HAVE_EXODUS_API

    // We write the file in the ExodusII format.
    ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

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

    // VariableGroup optimization means that DoFs iterate over each of
    // the 3 variables on each node before proceeding to the next
    // node.
    //                                     u_x, u_y, v
    const std::vector<Real> gold_vector = {-1,  3,   10,
                                           0,   1,   12,
                                           2,   0,   14,
                                           0,   1.5, 11,
                                           0.5, 1,   13,
                                           0.5, 1.5, 12,
                                           3,   4,   16,
                                           3,   1.5, 15,
                                           0.5, 4,   13};

    // Translation from node id to dof indexing order as gets done in
    // serial
    const std::vector<dof_id_type>
      node_reordering({0, 3, 1, 8, 5, 4, 6, 7, 2});

    Real tol = 1e-12;
    NumericVector<Number> & sys2_soln(*sys2.solution);
    for (const auto & node : read_mesh.node_ptr_range())
      {
        if (node->processor_id() != mesh.processor_id())
          continue;

        // We write out real, imaginary, and magnitude for complex
        // numbers
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
        const unsigned int v2 = 3;
#else
        const unsigned int v2 = 1;
#endif
        CPPUNIT_ASSERT_EQUAL(node->n_vars(0), 3*v2);

        const dof_id_type gold_i_ux = node_reordering[node->id()] * 3;
        const dof_id_type gold_i_uy = gold_i_ux + 1;
        const dof_id_type gold_i_v  = gold_i_uy + 1;

        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,0,0))),
                                     gold_vector[gold_i_ux], tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,v2,0))),
                                     gold_vector[gold_i_uy], tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,2*v2,0))),
                                     gold_vector[gold_i_v], tol);

        // Let's check imaginary parts and magnitude for good measure
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,1,0))),
                                     0.0, tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,2,0))),
                                     std::abs(gold_vector[gold_i_ux]), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,4,0))),
                                     0.0, tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,5,0))),
                                     std::abs(gold_vector[gold_i_uy]), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,7,0))),
                                     0.0, tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys2_soln(node->dof_number(0,8,0))),
                                     std::abs(gold_vector[gold_i_v]), tol);
#endif
      }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteVecAndScalar);
