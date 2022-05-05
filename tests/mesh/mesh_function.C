#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/dof_map.h>
#include <libmesh/system.h>
#include <libmesh/mesh_function.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/elem.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


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

Number trilinear_function (const Point & p,
                           const Parameters &,
                           const std::string &,
                           const std::string &)
{
  return 8*p(0) + 80*p(1) + 800*p(2);
}


class MeshFunctionTest : public CppUnit::TestCase
{
  /**
   * Tests for general MeshFunction capability.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshFunctionTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( test_subdomain_id_sets );
#endif
#if LIBMESH_DIM > 2
#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( test_p_level );
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

public:
  void setUp() {}

  void tearDown() {}

  // test that mesh function works correctly with subdomain id sets.
  void test_subdomain_id_sets()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_square (mesh,
                                         4, 4,
                                         0., 1.,
                                         0., 1.,
                                         QUAD4);

    // Set a subdomain id for all elements, based on location.
    for (auto & elem : mesh.active_element_ptr_range())
      {
        Point c = elem->vertex_average();
        elem->subdomain_id() =
          subdomain_id_type(c(0)*4) + subdomain_id_type(c(1)*4)*10;
      }

    // Add a discontinuous variable so we can easily see what side of
    // an interface we're querying
    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    unsigned int u_var = sys.add_variable("u", CONSTANT, MONOMIAL);

    es.init();
    sys.project_solution(trilinear_function, nullptr, es.parameters);

    MeshFunction mesh_function (sys.get_equation_systems(),
                                *sys.current_local_solution,
                                sys.get_dof_map(),
                                u_var);

    // Checkerboard pattern
    const std::set<subdomain_id_type> sbdids1 {0,2,11,13,20,22,31,33};
    const std::set<subdomain_id_type> sbdids2 {1,3,10,12,21,23,30,32};
    mesh_function.init();
    mesh_function.enable_out_of_mesh_mode(DenseVector<Number>());
    mesh_function.set_subdomain_ids(&sbdids1);

    DenseVector<Number> vec_values;
    const std::string dummy;

    // Make sure the MeshFunction's values interpolate the projected solution
    // at the nodes
    for (auto & elem : mesh.active_local_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        const Real expected_value =
          libmesh_real(trilinear_function(c, es.parameters, dummy, dummy));
        const std::vector<Point> offsets
          {{0,-1/8.}, {1/8.,0}, {0,1/8.}, {-1/8.,0}};
        for (Point offset : offsets)
          {
            const Point p = c + offset;
            mesh_function(p, 0, vec_values, &sbdids1);
            const Number retval1 = vec_values.empty() ? -12345 : vec_values(0);
            mesh_function(p, 0, vec_values, &sbdids2);
            const Number retval2 = vec_values.empty() ? -12345 : vec_values(0);
            mesh_function(c, 0, vec_values, nullptr);
            const Number retval3 = vec_values.empty() ? -12345 : vec_values(0);

            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(retval3), expected_value,
                                    TOLERANCE * TOLERANCE);

            if (sbdids1.count(elem->subdomain_id()))
              {
                CPPUNIT_ASSERT(!sbdids2.count(elem->subdomain_id()));
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(retval1), expected_value,
                                        TOLERANCE * TOLERANCE);

                mesh_function(c, 0, vec_values, &sbdids2);
                CPPUNIT_ASSERT(vec_values.empty());
              }
            else
              {
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(retval2), expected_value,
                                        TOLERANCE * TOLERANCE);

                mesh_function(c, 0, vec_values, &sbdids1);
                CPPUNIT_ASSERT(vec_values.empty());
              }
          }
      }
  }

  // test that mesh function works correctly with non-zero
  // Elem::p_level() values.
#ifdef LIBMESH_ENABLE_AMR
  void test_p_level()
  {
    LOG_UNIT_TEST;

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
    // std::make_unique doesn't like if we try to use {u_var} in-place?
    std::vector<unsigned int> variables {u_var};

    auto mesh_function =
      std::make_unique<MeshFunction>(sys.get_equation_systems(),
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

        LIBMESH_ASSERT_FP_EQUAL
          (libmesh_real(vec_values(0)),
           libmesh_real(mesh_function_value),
           TOLERANCE*TOLERANCE);
      }
  }
#endif // LIBMESH_ENABLE_AMR
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshFunctionTest );
