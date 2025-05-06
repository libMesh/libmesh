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
#include "libmesh/enum_xdr_mode.h"


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
#ifdef LIBMESH_HAVE_PETSC
  CPPUNIT_TEST( vectorMeshFunctionLagrange );
  CPPUNIT_TEST( vectorMeshFunctionNedelec );
  CPPUNIT_TEST( vectorMeshFunctionRaviartThomas );
#endif // LIBMESH_HAVE_PETSC
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

            LIBMESH_ASSERT_NUMBERS_EQUAL
              (retval3, expected_value, TOLERANCE * TOLERANCE);

            if (sbdids1.count(elem->subdomain_id()))
              {
                CPPUNIT_ASSERT(!sbdids2.count(elem->subdomain_id()));
                LIBMESH_ASSERT_NUMBERS_EQUAL
                  (retval1, expected_value, TOLERANCE * TOLERANCE);

                mesh_function(c, 0, vec_values, &sbdids2);
                CPPUNIT_ASSERT(vec_values.empty());
              }
            else
              {
                LIBMESH_ASSERT_NUMBERS_EQUAL
                  (retval2, expected_value, TOLERANCE * TOLERANCE);

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

        LIBMESH_ASSERT_NUMBERS_EQUAL
          (vec_values(0), mesh_function_value, TOLERANCE*TOLERANCE);
      }
  }
#endif // LIBMESH_ENABLE_AMR

  // Tests the projection of Lagrange Vectors using MeshFunction
  void vectorMeshFunctionLagrange()
  {
    LOG_UNIT_TEST;

    // Reading mesh and solution infromation from XDA files
    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read("solutions/lagrange_vec_solution_mesh.xda");
    EquationSystems es(mesh);
    es.read("solutions/lagrange_vec_solution.xda", READ,
            EquationSystems::READ_HEADER |
                EquationSystems::READ_DATA |
                EquationSystems::READ_ADDITIONAL_DATA);
    es.update();

    // Pulling the correct system and variable infromation from
    // the XDA files (the default system name is "nl0")
    System & sys = es.get_system<System>("nl0");
    std::unique_ptr<NumericVector<Number>> mesh_function_vector =
        NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize(*mesh_function_vector);

    // Setting up libMesh::MeshFunction
    MeshFunction mesh_function(es, *mesh_function_vector,
                                sys.get_dof_map(), 0);
    mesh_function.init();

    // Defining input parameters for MeshFunction::operator()
    DenseVector<Gradient> output;
    const std::set<subdomain_id_type> * subdomain_ids = nullptr;
    const Point & p = Point(0.5, 0.5);

    // Suppling the Lagrange Vec value at center of mesh to output
    (mesh_function)(p, 0.0, output, subdomain_ids);

    // Expected value at center mesh
    Gradient output_expected = Gradient(0.100977281077292,0.201954562154583);

    LIBMESH_ASSERT_FP_EQUAL(output(0)(0), output_expected(0),
                            TOLERANCE * TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(output(0)(1), output_expected(1),
                            TOLERANCE * TOLERANCE);
  }

  // Tests the projection of Nedelec Vectors using MeshFunction
  void vectorMeshFunctionNedelec()
  {
    LOG_UNIT_TEST;

    // Reading mesh and solution infromation from XDA files
    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read("solutions/nedelec_one_solution_mesh.xda");
    EquationSystems es(mesh);
    es.read("solutions/nedelec_one_solution.xda", READ,
            EquationSystems::READ_HEADER |
                EquationSystems::READ_DATA |
                EquationSystems::READ_ADDITIONAL_DATA);
    es.update();

    // Pulling the correct system and variable infromation from
    // the XDA files (the default system name is "nl0")
    System & sys = es.get_system<System>("nl0");
    std::unique_ptr<NumericVector<Number>> mesh_function_vector =
        NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize(*mesh_function_vector);

    // Setting up libMesh::MeshFunction
    MeshFunction mesh_function(es, *mesh_function_vector,
                                sys.get_dof_map(), 0);
    mesh_function.init();

    // Defining input parameters for MeshFunction::operator()
    DenseVector<Gradient> output;
    const std::set<subdomain_id_type> * subdomain_ids = nullptr;
    const Point & p = Point(0.5, 0.5);

    // Suppling the Nedelec One value at center of mesh to output
    (mesh_function)(p, 0.0, output, subdomain_ids);

    // Expected value at center mesh
    Gradient output_expected = Gradient(0.0949202883998996,-0.0949202883918033);

    LIBMESH_ASSERT_FP_EQUAL(output(0)(0), output_expected(0),
                            TOLERANCE * TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(output(0)(1), output_expected(1),
                            TOLERANCE * TOLERANCE);
  }

  // Tests the projection of Raviart Thomas Vectors using MeshFunction
  void vectorMeshFunctionRaviartThomas()
  {
    LOG_UNIT_TEST;

    // Reading mesh and solution infromation from XDA files
    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read("solutions/raviart_thomas_solution_mesh.xda");
    EquationSystems es(mesh);
    es.read("solutions/raviart_thomas_solution.xda", READ,
            EquationSystems::READ_HEADER |
                EquationSystems::READ_DATA |
                EquationSystems::READ_ADDITIONAL_DATA);
    es.update();

    // Pulling the correct system and variable infromation from
    // the XDA files (the default system name is "nl0")
    System & sys = es.get_system<System>("nl0");
    std::unique_ptr<NumericVector<Number>> mesh_function_vector =
        NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize(*mesh_function_vector);

    // Setting up libMesh::MeshFunction
    MeshFunction mesh_function(es, *mesh_function_vector,
                                sys.get_dof_map(), 0);
    mesh_function.init();

    // Defining input parameters for MeshFunction::operator()
    DenseVector<Gradient> output;
    const std::set<subdomain_id_type> * subdomain_ids = nullptr;
    const Point & p = Point(0.5, 0.5);

    // Suppling the Raviart Thomas value at center of mesh to output
    (mesh_function)(p, 0.0, output, subdomain_ids);

    // Expected value at center mesh
    Gradient output_expected = Gradient(0.0772539939808116,-0.0772537479511396);

    LIBMESH_ASSERT_FP_EQUAL(output(0)(0), output_expected(0),
                            TOLERANCE * TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(output(0)(1), output_expected(1),
                            TOLERANCE * TOLERANCE);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(MeshFunctionTest);
