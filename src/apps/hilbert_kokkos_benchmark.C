// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "L2system.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/diff_solver.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_solver.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/steady_solver.h"
#include "libmesh/string_to_enum.h"

#include <chrono>
#include <iomanip>
#include <string>
#include <vector>

#ifdef LIBMESH_HAVE_KOKKOS
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__
#endif

using namespace libMesh;

namespace
{

enum class GoalKind
{
  analytic,
  fem
};

struct BenchmarkOptions
{
  unsigned int nx = 256;
  unsigned int ny = 256;
  unsigned int repeats = 1;
  unsigned int hilbert_order = 0;
  unsigned int fe_order = 1;
  unsigned int max_linear_iterations = 5000;
  std::string fe_family = "LAGRANGE";
  std::string calc_expression = "sin(pi*x) + 0.25*y";
  std::string input_expression = "sin(pi*x) + 0.5*y";
  ElemType elem_type = QUAD4;
  Real linear_tolerance = 1.e-10;
  Real fdm_eps = 1.e-7;
  SolverType solver_type = CG;
  PreconditionerType preconditioner_type = JACOBI_PRECOND;
  GoalKind goal_kind = GoalKind::analytic;
};

struct BenchmarkResult
{
  std::vector<Number> solution;
  Real average_solve_seconds = 0.;
  Real average_plan_seconds = 0.;
  Real average_assembly_seconds = 0.;
  Real average_solver_seconds = 0.;
  Real average_total_seconds = 0.;
  dof_id_type n_dofs = 0;
};

#ifdef LIBMESH_HAVE_KOKKOS
struct KokkosScope
{
  KokkosScope(int & argc, char ** & argv)
  {
    ::Kokkos::initialize(argc, argv);
  }

  ~KokkosScope()
  {
    ::Kokkos::finalize();
  }
};
#endif

void usage_error(const char * progname)
{
  libMesh::out
    << "Options: " << progname << '\n'
    << " --nx              n         number of elements in x         [default: 256]\n"
    << " --ny              n         number of elements in y         [default: 256]\n"
    << " --elem-type       type      element type                    [default: QUAD4]\n"
    << " --goal            kind      analytic|fem                    [default: analytic]\n"
    << " --calc            expr      projection goal expression\n"
    << "                                default analytic: sin(pi*x) + 0.25*y\n"
    << "                                default fem:      u*u + x - 0.25*y\n"
    << " --input-func      expr      input field for fem goal        [default: sin(pi*x) + 0.5*y]\n"
    << " --family          fam       output FE family                [default: LAGRANGE]\n"
    << " --order           p         output FE order                 [default: 1]\n"
    << " --hilbert         p         Hilbert order                   [default: 0]\n"
    << " --fdm_eps         eps       fallback finite-difference eps  [default: 1e-7]\n"
    << " --linear_tol      tol       linear solve tolerance          [default: 1e-10]\n"
    << " --max_its         n         linear max iterations           [default: 5000]\n"
    << " --solver          type      solver type                     [default: CG]\n"
    << " --pc              type      preconditioner type             [default: JACOBI_PRECOND]\n"
    << " --repeats         n         number of fresh runs to average [default: 1]\n"
    << std::endl;

  std::exit(1);
}

GoalKind parse_goal_kind(const std::string & goal_string)
{
  if (goal_string == "analytic")
    return GoalKind::analytic;
  if (goal_string == "fem")
    return GoalKind::fem;

  libmesh_error_msg("Unsupported --goal value '" << goal_string << "'. Use analytic or fem.");
}

void configure_hilbert_system(HilbertSystem & sys,
                              const BenchmarkOptions & options,
                              const bool use_kokkos)
{
  sys.hilbert_order() = options.hilbert_order;
  sys.fe_family() = options.fe_family;
  sys.fe_order() = options.fe_order;
  sys.use_kokkos_backend(use_kokkos);
  sys.set_fdm_eps(options.fdm_eps);
  sys.time_solver = std::make_unique<SteadySolver>(sys);
}

void configure_linear_solver(HilbertSystem & sys,
                             const BenchmarkOptions & options)
{
  DiffSolver & solver = *sys.time_solver->diff_solver();
  solver.quiet = true;
  solver.verbose = false;
  solver.relative_step_tolerance = 1.e-12;

  sys.parameters.set<unsigned int>("linear solver maximum iterations") =
    options.max_linear_iterations;
  sys.parameters.set<Real>("linear solver tolerance") =
    options.linear_tolerance;

  auto * linear_solver = sys.get_linear_solver();
  linear_solver->set_solver_type(options.solver_type);
  linear_solver->set_preconditioner_type(options.preconditioner_type);
}

BenchmarkResult solve_projection_once(const Parallel::Communicator & comm,
                                      const BenchmarkOptions & options,
                                      const bool use_kokkos)
{
  Mesh mesh(comm, 2);
  MeshTools::Generation::build_square(mesh,
                                      options.nx,
                                      options.ny,
                                      0.,
                                      1.,
                                      0.,
                                      1.,
                                      options.elem_type);

  EquationSystems es(mesh);
  ExplicitSystem * input_system = nullptr;
  if (options.goal_kind == GoalKind::fem)
    {
      input_system = &es.add_system<ExplicitSystem>("input");
      input_system->add_variable("u", FIRST, LAGRANGE);
    }

  HilbertSystem & projection = es.add_system<HilbertSystem>("projection");
  configure_hilbert_system(projection, options, use_kokkos);
  projection.input_system = input_system;
  es.init();

  if (options.goal_kind == GoalKind::analytic)
    {
      ParsedFunction<Number> goal(options.calc_expression);
      projection.set_goal_func(goal);
    }
  else
    {
      libmesh_assert(input_system);
      ParsedFunction<Number> input_function(options.input_expression);
      input_system->project_solution(&input_function);

      ParsedFEMFunction<Number> goal(*input_system, options.calc_expression);
      projection.set_goal_func(goal);
    }

  configure_linear_solver(projection, options);

  const auto start = std::chrono::steady_clock::now();
  projection.solve();
  const auto stop = std::chrono::steady_clock::now();

  BenchmarkResult result;
  result.average_solve_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(stop - start).count();
  result.average_total_seconds = result.average_solve_seconds;
  if (use_kokkos)
    {
      const auto & timing = projection.last_kokkos_timing();
      result.average_plan_seconds = timing.plan_seconds;
      result.average_assembly_seconds = timing.assembly_seconds;
      result.average_solver_seconds = timing.solve_seconds;
      result.average_total_seconds = timing.total_seconds;
    }
  result.n_dofs = projection.n_dofs();
  projection.solution->localize(result.solution);
  return result;
}

BenchmarkResult solve_projection(const Parallel::Communicator & comm,
                                 const BenchmarkOptions & options,
                                 const bool use_kokkos)
{
  BenchmarkResult result;

  for (unsigned int repeat = 0; repeat != options.repeats; ++repeat)
    {
      auto single = solve_projection_once(comm, options, use_kokkos);
      result.average_solve_seconds += single.average_solve_seconds;
      result.average_plan_seconds += single.average_plan_seconds;
      result.average_assembly_seconds += single.average_assembly_seconds;
      result.average_solver_seconds += single.average_solver_seconds;
      result.average_total_seconds += single.average_total_seconds;
      result.n_dofs = single.n_dofs;
      result.solution = std::move(single.solution);
    }

  result.average_solve_seconds /= options.repeats;
  result.average_plan_seconds /= options.repeats;
  result.average_assembly_seconds /= options.repeats;
  result.average_solver_seconds /= options.repeats;
  result.average_total_seconds /= options.repeats;
  return result;
}

void
print_kokkos_phase_diagnostics(const BenchmarkResult & result)
{
  if (result.average_total_seconds <= 0.)
    return;

  const Real plan_fraction = result.average_plan_seconds / result.average_total_seconds;
  const Real assembly_fraction = result.average_assembly_seconds / result.average_total_seconds;
  const Real solver_fraction = result.average_solver_seconds / result.average_total_seconds;
  const Real accounted_fraction = plan_fraction + assembly_fraction + solver_fraction;
  const Real other_fraction = std::max<Real>(0., 1. - accounted_fraction);

  const char * dominant_name = "plan";
  Real dominant_fraction = plan_fraction;
  if (assembly_fraction > dominant_fraction)
    {
      dominant_name = "assembly";
      dominant_fraction = assembly_fraction;
    }
  if (solver_fraction > dominant_fraction)
    {
      dominant_name = "solver";
      dominant_fraction = solver_fraction;
    }
  if (other_fraction > dominant_fraction)
    {
      dominant_name = "other";
      dominant_fraction = other_fraction;
    }

  libMesh::out << "Kokkos phase fractions:"
               << " plan=" << 100. * plan_fraction << '%'
               << " assembly=" << 100. * assembly_fraction << '%'
               << " solver=" << 100. * solver_fraction << '%'
               << " other=" << 100. * other_fraction << '%'
               << std::endl;
  libMesh::out << "Dominant Kokkos phase: "
               << dominant_name
               << " (" << 100. * dominant_fraction << "% of total)"
               << std::endl;
}

void compute_difference_metrics(const std::vector<Number> & host_solution,
                                const std::vector<Number> & kokkos_solution,
                                Real & max_abs_host,
                                Real & max_abs_diff)
{
  libmesh_assert_equal_to(host_solution.size(), kokkos_solution.size());

  max_abs_host = 0.;
  max_abs_diff = 0.;

  for (const auto i : index_range(host_solution))
    {
      max_abs_host = std::max(max_abs_host, std::abs(libmesh_real(host_solution[i])));
      max_abs_diff = std::max(max_abs_diff,
                              std::abs(libmesh_real(host_solution[i] - kokkos_solution[i])));
    }
}

BenchmarkOptions parse_options()
{
  BenchmarkOptions options;

  if (libMesh::on_command_line("--help"))
    usage_error("hilbert_kokkos_benchmark");

  options.nx = libMesh::command_line_next("--nx", options.nx);
  options.ny = libMesh::command_line_next("--ny", options.ny);
  options.repeats = libMesh::command_line_next("--repeats", options.repeats);
  options.hilbert_order = libMesh::command_line_next("--hilbert", options.hilbert_order);
  options.fe_order = libMesh::command_line_next("--order", options.fe_order);
  options.fe_family = libMesh::command_line_next("--family", options.fe_family);
  options.input_expression = libMesh::command_line_next("--input-func", options.input_expression);
  options.elem_type =
    Utility::string_to_enum<ElemType>(libMesh::command_line_next("--elem-type",
                                                                 std::string("QUAD4")));
  options.fdm_eps = libMesh::command_line_next("--fdm_eps", options.fdm_eps);
  options.linear_tolerance = libMesh::command_line_next("--linear_tol", options.linear_tolerance);
  options.max_linear_iterations =
    libMesh::command_line_next("--max_its", options.max_linear_iterations);
  options.solver_type =
    Utility::string_to_enum<SolverType>(libMesh::command_line_next("--solver",
                                                                   std::string("CG")));
  options.preconditioner_type =
    Utility::string_to_enum<PreconditionerType>(
      libMesh::command_line_next("--pc", std::string("JACOBI_PRECOND")));
  options.goal_kind =
    parse_goal_kind(libMesh::command_line_next("--goal", std::string("analytic")));

  const std::string default_calc =
    options.goal_kind == GoalKind::analytic ?
      std::string("sin(pi*x) + 0.25*y") :
      std::string("u*u + x - 0.25*y");
  options.calc_expression = libMesh::command_line_next("--calc", default_calc);

  libmesh_error_msg_if(options.nx == 0 || options.ny == 0,
                       "--nx and --ny must both be positive");
  libmesh_error_msg_if(options.repeats == 0,
                       "--repeats must be positive");

  return options;
}

} // namespace

int main(int argc, char ** argv)
{
#ifdef LIBMESH_HAVE_KOKKOS
  KokkosScope kokkos_scope(argc, argv);
#endif
  LibMeshInit init(argc, argv);

#ifndef LIBMESH_HAVE_KOKKOS
  libmesh_error_msg("hilbert_kokkos_benchmark requires a libMesh build with Kokkos enabled");
#endif
#ifndef LIBMESH_HAVE_PETSC
  libmesh_error_msg("hilbert_kokkos_benchmark requires a libMesh build with PETSc enabled");
#endif
#ifndef LIBMESH_HAVE_FPARSER
  libmesh_error_msg("hilbert_kokkos_benchmark requires a libMesh build with FPARSER enabled");
#endif
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  libmesh_error_msg("hilbert_kokkos_benchmark does not support complex Number builds");
#endif

  const auto options = parse_options();

  libMesh::out << std::setprecision(std::numeric_limits<Real>::max_digits10);
  libMesh::out << "Running Hilbert benchmark with"
               << " nx=" << options.nx
               << " ny=" << options.ny
               << " elem_type=" << Utility::enum_to_string(options.elem_type)
               << " goal=" << (options.goal_kind == GoalKind::analytic ? "analytic" : "fem")
               << " family=" << options.fe_family
               << " order=" << options.fe_order
               << " hilbert=" << options.hilbert_order
               << " solver=" << Utility::enum_to_string(options.solver_type)
               << " pc=" << Utility::enum_to_string(options.preconditioner_type)
               << " repeats=" << options.repeats
               << std::endl;

  libMesh::out << "Starting host projection" << std::endl;
  const auto host_result = solve_projection(init.comm(), options, false);

  libMesh::out << "Starting Kokkos projection" << std::endl;
  const auto kokkos_result = solve_projection(init.comm(), options, true);

  libmesh_assert_equal_to(host_result.n_dofs, kokkos_result.n_dofs);
  libmesh_assert_equal_to(host_result.solution.size(), kokkos_result.solution.size());

  Real max_abs_host = 0.;
  Real max_abs_diff = 0.;
  compute_difference_metrics(host_result.solution,
                             kokkos_result.solution,
                             max_abs_host,
                             max_abs_diff);

  libMesh::out << "Degrees of freedom: " << host_result.n_dofs << std::endl;
  libMesh::out << "Host solve time:    " << host_result.average_solve_seconds << " s" << std::endl;
  libMesh::out << "Kokkos solve time:  " << kokkos_result.average_solve_seconds << " s" << std::endl;
  libMesh::out << "Kokkos plan time:   " << kokkos_result.average_plan_seconds << " s" << std::endl;
  libMesh::out << "Kokkos assembly:    " << kokkos_result.average_assembly_seconds << " s" << std::endl;
  libMesh::out << "Kokkos solver:      " << kokkos_result.average_solver_seconds << " s" << std::endl;
  libMesh::out << "Kokkos total time:  " << kokkos_result.average_total_seconds << " s" << std::endl;
  print_kokkos_phase_diagnostics(kokkos_result);

  if (kokkos_result.average_solve_seconds > 0.)
    libMesh::out << "Host/Kokkos ratio:  "
                 << host_result.average_solve_seconds / kokkos_result.average_solve_seconds
                 << std::endl;

  libMesh::out << "Max |host|:         " << max_abs_host << std::endl;
  libMesh::out << "Max |host-kokkos|:  " << max_abs_diff << std::endl;

  return 0;
}
