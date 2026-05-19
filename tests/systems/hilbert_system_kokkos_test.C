#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <libmesh/equation_systems.h>
#include <libmesh/enum_preconditioner_type.h>
#include <libmesh/enum_solver_type.h>
#include <libmesh/explicit_system.h>
#include <libmesh/linear_solver.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/parsed_fem_function.h>
#include <libmesh/parsed_function.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/steady_solver.h>
#include <libmesh/diff_solver.h>

#include <algorithm>
#include <chrono>
#include <vector>

#include "../../src/apps/L2system.C"

using namespace libMesh;

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && \
    defined(LIBMESH_HAVE_FPARSER) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
namespace
{

constexpr Real projection_tolerance = 5.e-10;

struct TimedProjectionResult
{
  std::vector<Number> solution;
  Real elapsed_seconds = 0.;
};

void
configure_hilbert_system(HilbertSystem & sys, const bool use_kokkos)
{
  sys.hilbert_order() = 1;
  sys.fe_family() = "LAGRANGE";
  sys.fe_order() = 1;
  sys.use_kokkos_backend(use_kokkos);
  sys.use_exact_parsed_fem_host_path(true);
  sys.time_solver = std::make_unique<SteadySolver>(sys);
}

void
configure_linear_solver(HilbertSystem & sys)
{
  DiffSolver & solver = *sys.time_solver->diff_solver();
  solver.quiet = true;
  solver.verbose = false;
  solver.relative_step_tolerance = 1.e-12;

  sys.parameters.set<unsigned int>("linear solver maximum iterations") = 500;
  sys.parameters.set<Real>("linear solver tolerance") = 1.e-14;

  auto * linear_solver = sys.get_linear_solver();
  linear_solver->set_solver_type(CG);
  linear_solver->set_preconditioner_type(IDENTITY_PRECOND);
}

std::vector<Number>
localize_solution(const System & sys)
{
  std::vector<Number> values;
  sys.solution->localize(values);
  return values;
}

void
assert_solutions_close(const std::vector<Number> & host_solution,
                       const std::vector<Number> & kokkos_solution)
{
  CPPUNIT_ASSERT_EQUAL(host_solution.size(), kokkos_solution.size());

  Real max_abs_host = 0;
  Real max_abs_diff = 0;

  for (const auto i : index_range(host_solution))
    {
      max_abs_host = std::max(max_abs_host, std::abs(libmesh_real(host_solution[i])));
      max_abs_diff = std::max(max_abs_diff,
                              std::abs(libmesh_real(host_solution[i] - kokkos_solution[i])));
    }

  const Real scaled_tol = projection_tolerance * std::max<Real>(1., max_abs_host);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., max_abs_diff, scaled_tol);
}

template <typename SolveFunctor>
TimedProjectionResult
time_projection_solve(SolveFunctor && solve)
{
  const auto start = std::chrono::steady_clock::now();
  auto solution = solve();
  const auto stop = std::chrono::steady_clock::now();

  TimedProjectionResult result;
  result.solution = std::move(solution);
  result.elapsed_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(stop - start).count();
  return result;
}

void
report_projection_timing(const std::string & label,
                         const TimedProjectionResult & host_result,
                         const TimedProjectionResult & kokkos_result)
{
  libMesh::out << label
               << " host_time=" << host_result.elapsed_seconds << " s"
               << " kokkos_time=" << kokkos_result.elapsed_seconds << " s";

  if (kokkos_result.elapsed_seconds > 0.)
    libMesh::out << " host_over_kokkos="
                 << host_result.elapsed_seconds / kokkos_result.elapsed_seconds;

  libMesh::out << std::endl;
}

std::vector<Number>
solve_analytic_projection_impl(const bool use_kokkos)
{
  ReplicatedMesh mesh(*TestCommWorld);
  MeshTools::Generation::build_square(mesh,
                                      3,
                                      2,
                                      0.,
                                      1.,
                                      0.,
                                      1.,
                                      QUAD4);

  EquationSystems es(mesh);
  HilbertSystem & sys = es.add_system<HilbertSystem>("projection");
  configure_hilbert_system(sys, use_kokkos);
  es.init();

  ParsedFunction<Number> goal("sin(pi*x) + 0.25*y");
  sys.set_goal_func(goal);
  sys.set_fdm_eps(1.e-7);
  configure_linear_solver(sys);
  sys.solve();

  return localize_solution(sys);
}

TimedProjectionResult
solve_analytic_projection(const bool use_kokkos)
{
  return time_projection_solve([&]() { return solve_analytic_projection_impl(use_kokkos); });
}

std::vector<Number>
solve_parsed_fem_projection_impl(const bool use_kokkos)
{
  ReplicatedMesh mesh(*TestCommWorld);
  MeshTools::Generation::build_square(mesh,
                                      3,
                                      2,
                                      0.,
                                      1.,
                                      0.,
                                      1.,
                                      QUAD4);

  EquationSystems es(mesh);
  ExplicitSystem & input = es.add_system<ExplicitSystem>("input");
  input.add_variable("u", FIRST, LAGRANGE);

  HilbertSystem & sys = es.add_system<HilbertSystem>("projection");
  configure_hilbert_system(sys, use_kokkos);
  sys.input_system = &input;
  es.init();

  ParsedFunction<Number> input_projection("sin(pi*x) + 0.5*y");
  input.project_solution(&input_projection);

  ParsedFEMFunction<Number> goal(input, "u*u + x - 0.25*y");
  sys.set_goal_func(goal);
  sys.set_fdm_eps(1.e-7);
  configure_linear_solver(sys);
  sys.solve();

  return localize_solution(sys);
}

TimedProjectionResult
solve_parsed_fem_projection(const bool use_kokkos)
{
  return time_projection_solve([&]() { return solve_parsed_fem_projection_impl(use_kokkos); });
}

void
report_single_projection_timing(const std::string & label,
                                const TimedProjectionResult & result)
{
  libMesh::out << label << " time=" << result.elapsed_seconds << " s" << std::endl;
}

} // namespace
#endif

class HilbertSystemKokkosTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(HilbertSystemKokkosTest);
  CPPUNIT_TEST(testAnalyticParsedFunctionEquivalence);
  CPPUNIT_TEST(testParsedFEMFunctionEquivalence);
  CPPUNIT_TEST_SUITE_END();

  void testAnalyticParsedFunctionEquivalence()
  {
    LOG_UNIT_TEST;

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && \
    defined(LIBMESH_HAVE_FPARSER) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
    libMesh::out << "Starting analytic host solve" << std::endl;
    const auto host_solution = solve_analytic_projection(false);
    report_single_projection_timing("Hilbert analytic host projection", host_solution);
    libMesh::out << "Starting analytic kokkos solve" << std::endl;
    const auto kokkos_solution = solve_analytic_projection(true);
    report_single_projection_timing("Hilbert analytic kokkos projection", kokkos_solution);
    report_projection_timing("Hilbert analytic projection",
                             host_solution,
                             kokkos_solution);
    assert_solutions_close(host_solution.solution, kokkos_solution.solution);
#endif
  }

  void testParsedFEMFunctionEquivalence()
  {
    LOG_UNIT_TEST;

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && \
    defined(LIBMESH_HAVE_FPARSER) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
    libMesh::out << "Starting parsed FEM host solve" << std::endl;
    const auto host_solution = solve_parsed_fem_projection(false);
    report_single_projection_timing("Hilbert parsed FEM host projection", host_solution);
    libMesh::out << "Starting parsed FEM kokkos solve" << std::endl;
    const auto kokkos_solution = solve_parsed_fem_projection(true);
    report_single_projection_timing("Hilbert parsed FEM kokkos projection", kokkos_solution);
    report_projection_timing("Hilbert parsed FEM projection",
                             host_solution,
                             kokkos_solution);
    assert_solutions_close(host_solution.solution, kokkos_solution.solution);
#endif
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(HilbertSystemKokkosTest);
