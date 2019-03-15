// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/dof_map.h>
#include <libmesh/system.h>
#include <libmesh/mesh_function.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/elem.h>

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

Number projection_function (const Point & p,
                            const Parameters &,
                            const std::string &,
                            const std::string &)
{
  return
    cos(.5*libMesh::pi*p(0)) *
    sin(.5*libMesh::pi*p(1)) *
    cos(.5*libMesh::pi*p(2));
}

class MeshFunctionTest : public CppUnit::TestCase
{
  /**
   * Tests for general MeshFunction capability.
   */
public:
  CPPUNIT_TEST_SUITE( MeshFunctionTest );

#if LIBMESH_DIM > 1
#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( test_p_level );
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

public:
  void setUp() {}

  void tearDown() {}

  // test that mesh function works correctly with non-zero
  // Elem::p_level() values.
#ifdef LIBMESH_ENABLE_AMR
  void test_p_level()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                       5, 5, 5,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       HEX20);

    // Bump the p-level for all elements. We will create a system with
    // a FIRST, LAGRANGE variable and then use the additional p-level
    // to solve with quadratic elements. Note: set_p_level(1) actually
    // _increases_ the existing variable order by 1, it does not set
    // it to 1.
    for (auto & elem : mesh.active_element_ptr_range())
      elem->set_p_level(1);

    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    unsigned int u_var = sys.add_variable("u", FIRST, LAGRANGE);

    es.init();
    sys.project_solution(projection_function, nullptr, es.parameters);


    std::unique_ptr<NumericVector<Number>> mesh_function_vector
      = NumericVector<Number>::build(sys.comm());
    mesh_function_vector->init(sys.n_dofs(), sys.n_local_dofs(),
                               sys.get_dof_map().get_send_list(), false,
                               GHOSTED);

    sys.solution->localize(*mesh_function_vector,
                           sys.get_dof_map().get_send_list());

    // So the MeshFunction knows which variables to compute values for.
    std::vector<unsigned int> variables(1);
    variables[0] = u_var;

    auto mesh_function =
      libmesh_make_unique<MeshFunction>(sys.get_equation_systems(),
                                        *mesh_function_vector,
                                        sys.get_dof_map(),
                                        variables);

    mesh_function->init();
    mesh_function->set_point_locator_tolerance(0.0001);
    DenseVector<Number> vec_values;
    std::string dummy;

    // Make sure the MeshFunction's values interpolate the projected solution
    // at the nodes
    for (const auto & node : mesh.local_node_ptr_range())
      {
        (*mesh_function)(*node, /*time=*/ 0., vec_values);
        Number mesh_function_value =
          projection_function(*node,
                              es.parameters,
                              dummy,
                              dummy);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(vec_values(0)),
           libmesh_real(mesh_function_value),
           TOLERANCE*TOLERANCE);
      }
  }
#endif // LIBMESH_ENABLE_AMR
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshFunctionTest );
