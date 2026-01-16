#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/dof_map.h>
#include <libmesh/system.h>
#include <libmesh/elem.h>
#include <libmesh/numeric_vector.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>
#include <vector>

#include <libmesh/exodusII_io.h>

using namespace libMesh;

namespace
{

  Number baseline_linear_function(const Point &p,
                                  const Parameters &,
                                  const std::string &,
                                  const std::string &)
    {
      return 0.5 * p(0) + 0.25 * p(1);
    }

  Number circle_projection_function(const Point &p,
                                    const Parameters &,
                                    const std::string &,
                                    const std::string &)
    {
      return std::cos(0.5 * libMesh::pi * p(0)) *
            std::sin(0.5 * libMesh::pi * p(1)) *
            std::cos(0.5 * libMesh::pi * p(2));
    }

}

class ProjectSolutionTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(ProjectSolutionTest);
  CPPUNIT_TEST(test_partial_project_solution);
  CPPUNIT_TEST_SUITE_END();

private:
  void test_partial_project_solution()
  {
    LOG_UNIT_TEST;

#if LIBMESH_DIM < 2
    return;
#else
    ReplicatedMesh mesh(*TestCommWorld, /*dim=*/2);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/5, /*ny=*/5,
                                        /*xmin=*/-1., /*xmax=*/1.,
                                        /*ymin=*/-1., /*ymax=*/1.,
                                        QUAD4);

    EquationSystems es(mesh);
    System &sys = es.add_system<System>("ProjSys");

    // one dof per element
    const unsigned int u_var = sys.add_variable("u", CONSTANT, MONOMIAL);
    const unsigned int v_var = sys.add_variable("v", CONSTANT, MONOMIAL);

    es.init();

    // baseline
    sys.project_solution(baseline_linear_function, /*gptr=*/nullptr, es.parameters);
    // ExodusII_IO(mesh).write_equation_systems("before_project.e", es);

    const DofMap &dof_map = sys.get_dof_map();

    // collect elements whose vertex_average lies inside the circle
    const Real r2 = 0.5 * 0.5;

    std::vector<const Elem *> selected_elems;
    for (const auto &e : mesh.active_local_element_ptr_range())
      {
        const Point c = e->vertex_average();
        if (c(0) * c(0) + c(1) * c(1) <= r2)
          selected_elems.push_back(e);
      }

    // ConstElemRange expects mesh element iterators OR a vec_type*
    // Here we use the vector-backed constructor.
    ConstElemRange circle_range(&selected_elems);

    // project only on selected elements and only on u variable
    std::vector<unsigned int> vars_to_project = {u_var};
    sys.project_solution(circle_projection_function, /*gptr=*/nullptr, es.parameters, circle_range, vars_to_project);
    // ExodusII_IO(mesh).write_equation_systems("after_project.e", es);

    // Check that restricted projection only overwrites elements inside the circle.
    // Outside elements must keep the baseline projection.
    const Real tol = 1e-3;

    for (const auto &e : mesh.active_local_element_ptr_range())
      {
        const Point c = e->vertex_average();
        const bool inside = (c(0) * c(0) + c(1) * c(1) <= r2);

        std::vector<dof_id_type> u_indices, v_indices;
        dof_map.dof_indices(e, u_indices, u_var);
        dof_map.dof_indices(e, v_indices, v_var);

        const Real u_val = (*sys.solution)(u_indices[0]);
        const Real v_val = (*sys.solution)(v_indices[0]);

        const Real val_baseline = baseline_linear_function(c, es.parameters, "", "");
        const Real val_projected = circle_projection_function(c, es.parameters, "", "");
        if (inside)
          // Inside the circle, values should be equal to the projected value.
          LIBMESH_ASSERT_FP_EQUAL(u_val, val_projected, tol);
        else
          // Outside the circle, values must remain unchanged.
          LIBMESH_ASSERT_FP_EQUAL(u_val, val_baseline, tol);

        // v variable was never projected, should always equal baseline
        LIBMESH_ASSERT_FP_EQUAL(v_val, val_baseline, tol);
      }
#endif
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ProjectSolutionTest);
