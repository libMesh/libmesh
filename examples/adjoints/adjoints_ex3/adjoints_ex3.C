// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// <h1>Adjoints Example 3 - Stokes flow and Convection-Diffusion</h1>
// \author Vikram Garg
// \date 2012
//
// We solve a coupled Stokes + Convection Diffusion system in an H channel geometry
// with 2 inlets and 2 outlets. The QoI is the species
// flux from left outlet

// Channel Geometry:
//                              Wall
//  ----------------------------------------------------------------
// 0                                                                1
//  ----------------------------     -------------------------------
//                             |     |
//                        Wall |     | Wall
//                             |     |
// -----------------------------     -------------------------------
// 2                                                               2
// -----------------------------------------------------------------
//                              Wall

// The equations governing this flow are:
// Stokes: -VectorLaplacian(velocity) + grad(pressure) = vector(0)
// Convection-Diffusion: - dot(velocity, grad(concentration) ) + Laplacian(concentration) = 0

// The boundary conditions are:
// u_1(0) = -(y-2)*(y-3), u_2(0) = 0 ; u_1(1) = (y-2)*(y-3), u_2(1) = 0 ;
// u_1(walls) = 0, u_2(walls) = 0;
// C(0) = 1 ; C(1) = 0;
// grad(C) dot n (walls) = 0 ;
// grad(C) dot n (2) = 0 ; grad(C) dot n (3) = 0

// The QoI is:
// Q((u,v), C) = integral_{left_outlet}  - u * C ds

// The complete equal order adjoint QoI error estimate is: (Notation for derivatives: grad(C) = C,1 + C,2)
// Q(u) - Q(u_h) \leq
// |e(u_1)|_{H1} |e(u_1^*)|_{H1} +  |e(u_2)|_{H1} |e(u_2^*)|_{H1} +  (1/Pe) * |e(C)|_{H1} |e(C^*)|_{H1}  +
//  ||e(u_1,1)||_{L2} ||e(p^*)||_{L2} + ||e(u_2,2)||_{L2} ||e(p^*)||_{L2} + ||e(u_1,1^*)||_{L2} ||e(p)||_{L2} + ||e(u_2,2^*)||_{L2} ||e(p)||_{L2}  +
// ||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} + ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2}
// ||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} + ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2}
// = error_non_pressure + error_with_pressure + error_convection_diffusion_x + error_convection_diffusion_y

// C++ includes
#include <iostream>
#include <sys/time.h>
#include <iomanip>

// Libmesh includes
#include "libmesh/equation_systems.h"

#include "libmesh/twostep_time_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/steady_solver.h"

#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_diff_solver.h"

#include "libmesh/mesh.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/system.h"
#include "libmesh/system_norm.h"

#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/const_fem_function.h"
#include "libmesh/error_vector.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/qoi_set.h"
#include "libmesh/weighted_patch_recovery_error_estimator.h"

// Local includes
#include "coupled_system.h"
#include "domain.h"
#include "initial.h"
#include "femparameters.h"
#include "mysystems.h"
#include "output.h"
#include "H-qoi.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// Number output files

std::string numbered_filename(unsigned int t_step, // The timestep count
                              unsigned int a_step, // The adaptive step count
                              std::string solution_type, // primal or adjoint solve
                              std::string type,
                              std::string extension,
                              FEMParameters & param)
{
  std::ostringstream file_name;
  file_name << solution_type
            << ".out."
            << type
            << '.'
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step
            << '.'
            << std::setw(2)
            << std::setfill('0')
            << std::right
            << a_step
            << '.'
            << extension;

  if (param.output_bz2)
    file_name << ".bz2";
  else if (param.output_gz)
    file_name << ".gz";
  return file_name.str();
}

// Write tecplot, gmv and xda/xdr output

void write_output(EquationSystems & es,
                  unsigned int t_step, // The timestep count
                  unsigned int a_step, // The adaptive step count
                  std::string solution_type, // primal or adjoint solve
                  FEMParameters & param)
{
  MeshBase & mesh = es.get_mesh();

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      std::ostringstream file_name_gmv;
      file_name_gmv << solution_type
                    << ".out.gmv."
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step
                    << '.'
                    << std::setw(2)
                    << std::setfill('0')
                    << std::right
                    << a_step;

      GMVIO(mesh).write_equation_systems
        (file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.output_tecplot)
    {
      std::ostringstream file_name_tecplot;
      file_name_tecplot << solution_type
                        << ".out."
                        << std::setw(3)
                        << std::setfill('0')
                        << std::right
                        << t_step
                        << '.'
                        << std::setw(2)
                        << std::setfill('0')
                        << std::right
                        << a_step
                        << ".plt";

      TecplotIO(mesh).write_equation_systems
        (file_name_tecplot.str(), es);
    }
#endif

  if (param.output_xda || param.output_xdr)
    mesh.renumber_nodes_and_elements();
  if (param.output_xda)
    {
      mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param));
      es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xda", param),
               WRITE, EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
  if (param.output_xdr)
    {
      mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param));
      es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param),
               WRITE, EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
}

void write_output_headers(FEMParameters & param)
{
  // Only one processor needs to take care of headers.
  if (libMesh::global_processor_id() != 0)
    return;

  start_output(param.initial_timestep, "out_clocktime.m", "vector_clocktime");

  if (param.run_simulation)
    {
      start_output(param.initial_timestep, "out_activemesh.m", "table_activemesh");
      start_output(param.initial_timestep, "out_solvesteps.m", "table_solvesteps");

      if (param.timesolver_tolerance)
        {
          start_output(param.initial_timestep, "out_time.m", "vector_time");
          start_output(param.initial_timestep, "out_timesteps.m", "vector_timesteps");
        }
      if (param.steadystate_tolerance)
        start_output(param.initial_timestep, "out_changerate.m", "vector_changerate");
    }
}

void write_output_solvedata(EquationSystems & es,
                            unsigned int a_step,
                            unsigned int newton_steps,
                            unsigned int krylov_steps,
                            unsigned int tv_sec,
                            unsigned int tv_usec)
{
  MeshBase & mesh = es.get_mesh();
  unsigned int n_active_elem = mesh.n_active_elem();
  unsigned int n_active_dofs = es.n_active_dofs();

  if (mesh.processor_id() == 0)
    {
      // Write out the number of elements/dofs used
      std::ofstream activemesh ("out_activemesh.m",
                                std::ios_base::app | std::ios_base::out);
      activemesh.precision(17);
      activemesh << (a_step + 1) << ' '
                 << n_active_elem << ' '
                 << n_active_dofs << std::endl;

      // Write out the number of solver steps used
      std::ofstream solvesteps ("out_solvesteps.m",
                                std::ios_base::app | std::ios_base::out);
      solvesteps.precision(17);
      solvesteps << newton_steps << ' '
                 << krylov_steps << std::endl;

      // Write out the clock time
      std::ofstream clocktime ("out_clocktime.m",
                               std::ios_base::app | std::ios_base::out);
      clocktime.precision(17);
      clocktime << tv_sec << '.' << tv_usec << std::endl;
    }
}

void write_output_footers(FEMParameters & param)
{
  // Write footers on output .m files
  if (libMesh::global_processor_id() == 0)
    {
      std::ofstream clocktime ("out_clocktime.m",
                               std::ios_base::app | std::ios_base::out);
      clocktime << "];" << std::endl;

      if (param.run_simulation)
        {
          std::ofstream activemesh ("out_activemesh.m",
                                    std::ios_base::app | std::ios_base::out);
          activemesh << "];" << std::endl;

          std::ofstream solvesteps ("out_solvesteps.m",
                                    std::ios_base::app | std::ios_base::out);
          solvesteps << "];" << std::endl;

          if (param.timesolver_tolerance)
            {
              std::ofstream times ("out_time.m",
                                   std::ios_base::app | std::ios_base::out);
              times << "];" << std::endl;
              std::ofstream timesteps ("out_timesteps.m",
                                       std::ios_base::app | std::ios_base::out);
              timesteps << "];" << std::endl;
            }
          if (param.steadystate_tolerance)
            {
              std::ofstream changerate ("out_changerate.m",
                                        std::ios_base::app | std::ios_base::out);
              changerate << "];" << std::endl;
            }
        }

      // We're done, so write out a file (for e.g. Make to check)
      std::ofstream complete ("complete");
      complete << "complete" << std::endl;
    }
}

#if defined(LIBMESH_HAVE_GMV) || defined(LIBMESH_HAVE_TECPLOT_API)
void write_error(EquationSystems & es,
                 ErrorVector & error,
                 unsigned int t_number,
                 unsigned int a_number,
                 FEMParameters & param,
                 std::string error_type)
#else
  void write_error(EquationSystems &,
                   ErrorVector &,
                   unsigned int,
                   unsigned int,
                   FEMParameters &,
                   std::string)
#endif
{
#ifdef LIBMESH_HAVE_GMV
  if (param.write_gmv_error)
    {
      std::ostringstream error_gmv;
      error_gmv << "error.gmv."
                << std::setw(3)
                << std::setfill('0')
                << std::right
                << a_number
                << "."
                << std::setw(2)
                << std::setfill('0')
                << std::right
                << t_number
                << error_type;

      error.plot_error(error_gmv.str(), es.get_mesh());
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.write_tecplot_error)
    {
      std::ostringstream error_tecplot;
      error_tecplot << "error.plt."
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << a_number
                    << "."
                    << std::setw(2)
                    << std::setfill('0')
                    << std::right
                    << t_number
                    << error_type;

      error.plot_error(error_tecplot.str(), es.get_mesh());
    }
#endif
}

void read_output(EquationSystems & es,
                 unsigned int t_step,
                 unsigned int a_step,
                 std::string solution_type,
                 FEMParameters & param)
{
  MeshBase & mesh = es.get_mesh();

  std::string file_name_mesh, file_name_soln;
  // Look for ASCII files first
  if (param.output_xda)
    {
      file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param);
      file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xda", param);
    }
  else if (param.output_xdr)
    {
      file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param);
      file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param);
    }

  // Read in the mesh
  mesh.read(file_name_mesh);

  // And the stored solution
  es.read(file_name_soln, READ,
          EquationSystems::READ_HEADER |
          EquationSystems::READ_DATA |
          EquationSystems::READ_ADDITIONAL_DATA);

  // Put systems in a consistent state
  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).update();

  // Figure out the current time
  Real current_time = 0., current_timestep = 0.;

  if (param.timesolver_tolerance)
    {
      std::ifstream times ("out_time.m");
      std::ifstream timesteps ("out_timesteps.m");
      if (times.is_open() && timesteps.is_open())
        {
          // Read headers
          const unsigned int headersize = 25;
          char header[headersize];
          timesteps.getline (header, headersize);
          if (strcmp(header, "vector_timesteps = [") != 0)
            libmesh_error_msg("Bad header in out_timesteps.m:\n" << header);

          times.getline (header, headersize);
          if (strcmp(header, "vector_time = [") != 0)
            libmesh_error_msg("Bad header in out_time.m:\n" << header);

          // Read each timestep
          for (unsigned int i = 0; i != t_step; ++i)
            {
              if (!times.good())
                libmesh_error_msg("Error: File out_time.m is in non-good state.");
              times >> current_time;
              timesteps >> current_timestep;
            }
          // Remember to increment the last timestep; out_times.m
          // lists each *start* time
          current_time += current_timestep;
        }
      else
        libmesh_error_msg("Error opening out_time.m or out_timesteps.m");
    }
  else
    current_time = t_step * param.deltat;

  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).time = current_time;
}

void set_system_parameters(FEMSystem & system,
                           FEMParameters & param)
{
  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;

  // More desperate debugging options
  system.print_solution_norms = param.print_solution_norms;
  system.print_solutions      = param.print_solutions;
  system.print_residual_norms = param.print_residual_norms;
  system.print_residuals      = param.print_residuals;
  system.print_jacobian_norms = param.print_jacobian_norms;
  system.print_jacobians      = param.print_jacobians;

  system.print_element_solutions = param.print_element_solutions;
  system.print_element_residuals = param.print_element_residuals;

  // Solve this as a time-dependent or steady system
  if (param.transient)
    {
      UnsteadySolver *innersolver;
      if (param.timesolver_core == "euler2")
        {
          Euler2Solver *euler2solver =
            new Euler2Solver(system);

          euler2solver->theta = param.timesolver_theta;
          innersolver = euler2solver;
        }
      else if (param.timesolver_core == "euler")
        {
          EulerSolver *eulersolver =
            new EulerSolver(system);

          eulersolver->theta = param.timesolver_theta;
          innersolver = eulersolver;
        }
      else
        libmesh_error_msg("Don't recognize core TimeSolver type: " << param.timesolver_core);

      if (param.timesolver_tolerance)
        {
          TwostepTimeSolver *timesolver =
            new TwostepTimeSolver(system);

          timesolver->max_growth       = param.timesolver_maxgrowth;
          timesolver->target_tolerance = param.timesolver_tolerance;
          timesolver->upper_tolerance  = param.timesolver_upper_tolerance;
          timesolver->component_norm   = SystemNorm(param.timesolver_norm);

          timesolver->core_time_solver.reset(innersolver);
          system.time_solver.reset(timesolver);
        }
      else
        system.time_solver.reset(innersolver);
    }
  else
    system.time_solver.reset(new SteadySolver(system));

  system.time_solver->reduce_deltat_on_diffsolver_failure =
    param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver().reset(solver);
#else
      libmesh_error_msg("This example requires libMesh to be compiled with PETSc support.");
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver().reset(solver);

      solver->quiet                       = param.solver_quiet;
      solver->verbose                     = !param.solver_quiet;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
      solver->absolute_residual_tolerance = param.absolute_residual_tolerance;
      solver->relative_residual_tolerance = param.relative_residual_tolerance;
      solver->require_residual_reduction  = param.require_residual_reduction;
      solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      if (system.time_solver->reduce_deltat_on_diffsolver_failure)
        {
          solver->continue_after_max_iterations = true;
          solver->continue_after_backtrack_failure = true;
        }

      // And the linear solver options
      solver->max_linear_iterations       = param.max_linear_iterations;
      solver->initial_linear_tolerance    = param.initial_linear_tolerance;
      solver->minimum_linear_tolerance    = param.minimum_linear_tolerance;
    }
}

#ifdef LIBMESH_ENABLE_AMR

UniquePtr<MeshRefinement> build_mesh_refinement(MeshBase & mesh,
                                                FEMParameters & param)
{
  MeshRefinement * mesh_refinement = new MeshRefinement(mesh);
  mesh_refinement->coarsen_by_parents() = true;
  mesh_refinement->absolute_global_tolerance() = param.global_tolerance;
  mesh_refinement->nelem_target()      = param.nelem_target;
  mesh_refinement->refine_fraction()   = param.refine_fraction;
  mesh_refinement->coarsen_fraction()  = param.coarsen_fraction;
  mesh_refinement->coarsen_threshold() = param.coarsen_threshold;

  return UniquePtr<MeshRefinement>(mesh_refinement);
}

#endif // LIBMESH_ENABLE_AMR

// This function builds the Kelly error indicator. This indicator can be used
// for comparisons of adjoint and non-adjoint based error indicators
UniquePtr<ErrorEstimator> build_error_estimator(FEMParameters & /* param */)
{
  return UniquePtr<ErrorEstimator>(new KellyErrorEstimator);
}

// Functions to build the adjoint based error indicators
// The error_non_pressure and error_pressure constributions are estimated using
// the build_error_estimator_component_wise function below
UniquePtr<ErrorEstimator>
build_error_estimator_component_wise (FEMParameters & param,
                                      std::vector<std::vector<Real> > & term_weights,
                                      std::vector<FEMNormType> & primal_error_norm_type,
                                      std::vector<FEMNormType> & dual_error_norm_type)
{
  AdjointResidualErrorEstimator * adjoint_residual_estimator = new AdjointResidualErrorEstimator;

  // Both the primal and dual weights are going to be estimated using the patch recovery error estimator
  PatchRecoveryErrorEstimator * p1 = new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->primal_error_estimator().reset(p1);

  PatchRecoveryErrorEstimator * p2 = new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->dual_error_estimator().reset(p2);

  // Set the boolean for specifying whether we are reusing patches while building the patch recovery estimates
  p1->set_patch_reuse(param.patch_reuse);
  p2->set_patch_reuse(param.patch_reuse);

  // Using the user filled error norm type vector, we pass the type of norm to be used for
  // the error in each variable, we can have different types of norms for the primal and
  // dual variables
  std::size_t size = primal_error_norm_type.size();

  libmesh_assert_equal_to (size, dual_error_norm_type.size());
  for (std::size_t i = 0; i != size; ++i)
    {
      adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(i, primal_error_norm_type[i]);
      adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(i, dual_error_norm_type[i]);
    }

  // Now we set the right weights for each term in the error estimate, using the user provided
  // term_weights matrix
  libmesh_assert_equal_to (size, term_weights.size());
  for (std::size_t i = 0; i != term_weights.size(); ++i)
    {
      libmesh_assert_equal_to (size, term_weights[i].size());
      adjoint_residual_estimator->error_norm.set_weight(i, term_weights[i][i]);
      for (std::size_t j = 0; j != size; ++j)
        if (i != j)
          adjoint_residual_estimator->error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
    }

  return UniquePtr<ErrorEstimator>(adjoint_residual_estimator);
}

// The error_convection_diffusion_x and error_convection_diffusion_y are the nonlinear contributions which
// are computed using the build_weighted_error_estimator_component_wise below
UniquePtr<ErrorEstimator>
build_weighted_error_estimator_component_wise (FEMParameters & param,
                                               std::vector<std::vector<Real> > & term_weights,
                                               std::vector<FEMNormType> & primal_error_norm_type,
                                               std::vector<FEMNormType> & dual_error_norm_type,
                                               std::vector<FEMFunctionBase<Number> *> coupled_system_weight_functions)
{
  AdjointResidualErrorEstimator * adjoint_residual_estimator = new AdjointResidualErrorEstimator;

  // Using the user filled error norm type vector, we pass the type of norm to be used for
  // the error in each variable, we can have different types of norms for the primal and
  // dual variables

  WeightedPatchRecoveryErrorEstimator * p1 = new WeightedPatchRecoveryErrorEstimator;
  adjoint_residual_estimator->primal_error_estimator().reset(p1);

  PatchRecoveryErrorEstimator * p2 = new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->dual_error_estimator().reset(p2);

  p1->set_patch_reuse(param.patch_reuse);
  p2->set_patch_reuse(param.patch_reuse);

  // This is the critical difference with the build_error_estimate_component_wise
  p1->weight_functions.clear();

  // We pass the pointers to the user specified weight functions to the patch recovery
  // error estimator objects declared above
  std::size_t size = primal_error_norm_type.size();

  libmesh_assert(coupled_system_weight_functions.size() == size);
  libmesh_assert(dual_error_norm_type.size() == size);

  for (unsigned int i = 0; i != size; ++i)
    {
      p1->weight_functions.push_back(coupled_system_weight_functions[i]);
      adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(i, primal_error_norm_type[i]);
      adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(i, dual_error_norm_type[i]);
    }

  // Now we set the right weights for each term in the error estimate, using the user provided
  // term_weights matrix
  libmesh_assert_equal_to (size, term_weights.size());
  for (std::size_t i = 0; i != term_weights.size(); ++i)
    {
      libmesh_assert_equal_to (size, term_weights[i].size());
      adjoint_residual_estimator->error_norm.set_weight(i, term_weights[i][i]);
      for (std::size_t j = 0; j != size; ++j)
        if (i != j)
          adjoint_residual_estimator->error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
    }

  return UniquePtr<ErrorEstimator>(adjoint_residual_estimator);
}

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // This doesn't converge with Eigen BICGSTAB or with Trilinos for some reason...
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  libMesh::out << "Started " << argv[0] << std::endl;

  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("general.in");
    if (!i)
      libmesh_error_msg('[' << init.comm().rank() << "] Can't find general.in; exiting early.");
  }

  GetPot infile("general.in");

  // Read in parameters from the input file
  FEMParameters param(init.comm());
  param.read(infile);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh with the given dimension, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm(), param.dimension);

  // And an object to refine it
  UniquePtr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  libMesh::out << "Building mesh" << std::endl;

  if (!param.initial_timestep && param.run_simulation)
    build_domain(mesh, param);

  libMesh::out << "Building system" << std::endl;

  FEMSystem & system = build_system(equation_systems, infile, param);

  set_system_parameters(system, param);

  libMesh::out << "Initializing systems" << std::endl;

  // Create a QoI object, we will need this to compute the QoI
  {
    CoupledSystemQoI qoi;
    // Our QoI is computed on the side
    qoi.assemble_qoi_sides = true;
    system.attach_qoi(&qoi);
  }

  equation_systems.init ();

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  libMesh::out << "Starting adaptive loop" << std::endl << std::endl;

  // Count the number of steps used
  unsigned int a_step = 0;

  // Get the linear solver object to set the preconditioner reuse flag
  LinearSolver<Number> * linear_solver = system.get_linear_solver();

  // Adaptively solve the timestep
  for (; a_step <= param.max_adaptivesteps; ++a_step)
    {
      // We have to ensure that we are refining to either an error
      // tolerance or to a target number of elements, and, not to a
      // non-zero tolerance and a target number of elements
      // simultaneously
      if (param.global_tolerance > 0 && param.nelem_target > 0)
        {
          libMesh::out<<"We cant refine to both a non-zero tolerance and a target number of elements, EXITING adaptive loop. "<<std::endl<<std::endl;
          break;
        }

      // We dont want to reuse the preconditioner for the forward solve
      linear_solver->reuse_preconditioner(false);

      libMesh::out << "Adaptive step " << a_step << std::endl;

      libMesh::out << "Solving the forward problem" << std::endl;

      libMesh::out << "We have "
                   << mesh.n_active_elem()
                   << " active elements and "
                   << equation_systems.n_active_dofs()
                   << " active dofs."
                   << std::endl
                   << std::endl;

      // Solve the forward system
      system.solve();

      // Write the output
      write_output(equation_systems, 0, a_step, "primal", param);

      // Compute the QoI
      system.assemble_qoi();

      // We just call a postprocess here to set the variable computed_QoI to the value computed by assemble_qoi
      system.postprocess();

      // Get the value of the computed_QoI variable of the CoupledSystem class
      Number QoI_0_computed = (dynamic_cast<CoupledSystem &>(system)).get_QoI_value();

      libMesh::out << "The boundary QoI is "
                   << std::setprecision(17)
                   << QoI_0_computed
                   << std::endl
                   << std::endl;

      // Hard coded assert to ensure that the actual numbers we are getting are what they should be
      if (a_step == param.max_adaptivesteps)
        libmesh_assert_less(std::abs(QoI_0_computed - 0.0833), 2.5e-4);

      // You dont need to compute error estimates and refine at the last
      // adaptive step, only before that
      if (a_step!=param.max_adaptivesteps)
        {
          NumericVector<Number> & primal_solution = (*system.solution);

          // Make sure we get the contributions to the adjoint RHS from the sides
          system.assemble_qoi_sides = true;

          // We are about to solve the adjoint system, but before we
          // do this we see the same preconditioner flag to reuse the
          // preconditioner from the forward solver
          linear_solver->reuse_preconditioner(param.reuse_preconditioner);

          // Solve the adjoint system
          libMesh::out<< "Solving the adjoint problem" <<std::endl;
          system.adjoint_solve();

          // Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
          system.set_adjoint_already_solved(true);

          // To plot the adjoint solution, we swap it with the primal solution
          // and use the write_output function
          NumericVector<Number> & dual_solution = system.get_adjoint_solution();
          primal_solution.swap(dual_solution);

          write_output(equation_systems, 0, a_step, "adjoint", param);

          // Swap back
          primal_solution.swap(dual_solution);

          // We need the values of the parameters Pe from the
          // system for the adjoint error estimate
          Real Pe = (dynamic_cast<CoupledSystem &>(system)).get_Pe();

          // The total error is the sum: error = error_non_pressure +
          // error_with_pressure + ...
          // error_estimator_convection_diffusion_x +
          // error_estimator_convection_diffusion_y
          ErrorVector error;

          // We first construct the non-pressure contributions
          ErrorVector error_non_pressure;

          // First we build the norm_type_vector_non_pressure vectors and
          // weights_matrix_non_pressure matrix for the non-pressure term
          // error contributions
          std::vector<FEMNormType> primal_norm_type_vector_non_pressure;
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          primal_norm_type_vector_non_pressure.push_back(L2);
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);

          std::vector<FEMNormType> dual_norm_type_vector_non_pressure;
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          dual_norm_type_vector_non_pressure.push_back(L2);
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);

          std::vector<std::vector<Real> >
            weights_matrix_non_pressure(system.n_vars(),
                                        std::vector<Real>(system.n_vars(), 0.0));
          weights_matrix_non_pressure[0][0] = 1.;
          weights_matrix_non_pressure[1][1] = 1.;
          weights_matrix_non_pressure[3][3] = 1./Pe;

          // We build the error estimator to estimate the contributions
          // to the QoI error from the non pressure term
          UniquePtr<ErrorEstimator> error_estimator_non_pressure =
            build_error_estimator_component_wise (param,
                                                  weights_matrix_non_pressure,
                                                  primal_norm_type_vector_non_pressure,
                                                  dual_norm_type_vector_non_pressure);

          // Estimate the contributions to the QoI error from the non
          // pressure terms
          error_estimator_non_pressure->estimate_error(system, error_non_pressure);

          // Plot the estimated error from the non_pressure terms
          write_error(equation_systems, error_non_pressure, 0, a_step, param, "_non_pressure");

          // Now for the pressure contributions
          ErrorVector error_with_pressure;

          // Next we build the norm_type_vector_with_pressure vectors and
          // weights_matrix_with_pressure matrix for the pressure term
          // error contributions
          std::vector<FEMNormType> primal_norm_type_vector_with_pressure;
          primal_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
          primal_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
          primal_norm_type_vector_with_pressure.push_back(L2);
          primal_norm_type_vector_with_pressure.push_back(L2);

          std::vector<FEMNormType> dual_norm_type_vector_with_pressure;
          dual_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
          dual_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
          dual_norm_type_vector_with_pressure.push_back(L2);
          dual_norm_type_vector_with_pressure.push_back(L2);

          std::vector<std::vector<Real> >
            weights_matrix_with_pressure (system.n_vars(),
                                          std::vector<Real>(system.n_vars(), 0.0));
          weights_matrix_with_pressure[0][2] = 1.;
          weights_matrix_with_pressure[1][2] = 1.;
          weights_matrix_with_pressure[2][0] = 1.;
          weights_matrix_with_pressure[2][1] = 1.;

          // We build the error estimator to estimate the contributions
          // to the QoI error from the pressure term
          UniquePtr<ErrorEstimator> error_estimator_with_pressure =
            build_error_estimator_component_wise (param,
                                                  weights_matrix_with_pressure,
                                                  primal_norm_type_vector_with_pressure,
                                                  dual_norm_type_vector_with_pressure);

          // Estimate the contributions to the QoI error from the pressure terms
          error_estimator_with_pressure->estimate_error(system, error_with_pressure);

          // Plot the error due to the pressure terms
          write_error(equation_systems, error_with_pressure, 0, a_step, param, "_with_pressure");

          // Now for the convection diffusion term errors (in the x and y directions)

          ErrorVector error_convection_diffusion_x;

          // The norm type vectors and weights matrix for the convection_diffusion_x errors
          std::vector<FEMNormType> primal_norm_type_vector_convection_diffusion_x;
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(H1_X_SEMINORM);

          std::vector<FEMNormType> dual_norm_type_vector_convection_diffusion_x;
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          // Note that we need the error of the dual concentration in L2
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);

          std::vector<std::vector<Real> >
            weights_matrix_convection_diffusion_x (system.n_vars(),
                                                   std::vector<Real>(system.n_vars(), 0.0));
          weights_matrix_convection_diffusion_x[0][3] = 1.;
          weights_matrix_convection_diffusion_x[3][3] = 1.;

          // We will also have to build and pass the weight functions to the weighted patch recovery estimators

          // We pass the identity function as weights to error entries that the above matrix will scale to 0.
          ConstFEMFunction<Number> identity(1);

          // Declare object of class CoupledFEMFunctionsx, the definition of the function contains the weight
          // to be applied to the relevant terms
          // For ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2} term, returns C,1_h
          CoupledFEMFunctionsx convdiffx0(system, 0);
          // For ||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} term, returns (u_1)_h
          CoupledFEMFunctionsx convdiffx3(system, 3);

          // Make a vector of pointers to these objects
          std::vector<FEMFunctionBase<Number> *> coupled_system_weight_functions_x;
          coupled_system_weight_functions_x.push_back(&convdiffx0);
          coupled_system_weight_functions_x.push_back(&identity);
          coupled_system_weight_functions_x.push_back(&identity);
          coupled_system_weight_functions_x.push_back(&convdiffx3);

          // Build the error estimator to estimate the contributions
          // to the QoI error from the convection diffusion x term
          UniquePtr<ErrorEstimator> error_estimator_convection_diffusion_x =
            build_weighted_error_estimator_component_wise (param,
                                                           weights_matrix_convection_diffusion_x,
                                                           primal_norm_type_vector_convection_diffusion_x,
                                                           dual_norm_type_vector_convection_diffusion_x,
                                                           coupled_system_weight_functions_x);

          // Estimate the contributions to the QoI error from the
          // convection diffusion x term
          error_estimator_convection_diffusion_x->estimate_error(system, error_convection_diffusion_x);

          // Plot this error
          write_error (equation_systems,
                       error_convection_diffusion_x,
                       0,
                       a_step,
                       param,
                       "_convection_diffusion_x");

          // Now for the y direction terms
          ErrorVector error_convection_diffusion_y;

          // The norm type vectors and weights matrix for the convection_diffusion_x errors
          std::vector<FEMNormType> primal_norm_type_vector_convection_diffusion_y;
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(H1_Y_SEMINORM);

          std::vector<FEMNormType> dual_norm_type_vector_convection_diffusion_y;
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          // Note that we need the error of the dual concentration in L2
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);

          std::vector<std::vector<Real> >
            weights_matrix_convection_diffusion_y (system.n_vars(),
                                                   std::vector<Real>(system.n_vars(), 0.0));
          weights_matrix_convection_diffusion_y[1][3] = 1.;
          weights_matrix_convection_diffusion_y[3][3] = 1.;

          // For ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2} term, returns C,2_h
          CoupledFEMFunctionsy convdiffy1(system, 1);
          // For ||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} term, returns (u_2)_h
          CoupledFEMFunctionsy convdiffy3(system, 3);

          // Make a vector of pointers to these objects
          std::vector<FEMFunctionBase<Number> *> coupled_system_weight_functions_y;
          coupled_system_weight_functions_y.push_back(&identity);
          coupled_system_weight_functions_y.push_back(&convdiffy1);
          coupled_system_weight_functions_y.push_back(&identity);
          coupled_system_weight_functions_y.push_back(&convdiffy3);

          // Build the error estimator to estimate the contributions
          // to the QoI error from the convection diffsion y term
          UniquePtr<ErrorEstimator> error_estimator_convection_diffusion_y =
            build_weighted_error_estimator_component_wise (param,
                                                           weights_matrix_convection_diffusion_y,
                                                           primal_norm_type_vector_convection_diffusion_y,
                                                           dual_norm_type_vector_convection_diffusion_y,
                                                           coupled_system_weight_functions_y);

          // Estimate the contributions to the QoI error from the
          // convection diffusion y terms
          error_estimator_convection_diffusion_y->estimate_error(system, error_convection_diffusion_y);

          // Plot this error
          write_error(equation_systems, error_convection_diffusion_y, 0, a_step, param, "_convection_diffusion_y");

          if (param.indicator_type == "adjoint_residual")
            {
              error.resize(error_non_pressure.size());

              // Now combine the contribs from the pressure and non
              // pressure terms to get the complete estimate
              for (std::size_t i = 0; i < error.size(); i++)
                error[i] = error_non_pressure[i] + error_with_pressure[i] + error_convection_diffusion_x[i] + error_convection_diffusion_y[i];
            }
          else
            {
              libMesh::out << "Using Kelly Estimator" << std::endl;

              // Build the Kelly error estimator
              UniquePtr<ErrorEstimator> error_estimator = build_error_estimator(param);

              // Estimate the error
              error_estimator->estimate_error(system, error);
            }

          write_error(equation_systems, error, 0, a_step, param, "_total");

          // We have to refine either based on reaching an error
          // tolerance or a number of elements target, which should be
          // verified above

          // Otherwise we flag elements by error tolerance or nelem
          // target

          if (param.refine_uniformly)
            {
              mesh_refinement->uniformly_refine(1);
            }
          else if (param.global_tolerance >= 0. && param.nelem_target == 0.) // Refine based on reaching an error tolerance
            {
              mesh_refinement->flag_elements_by_error_tolerance (error);

              // Carry out the adaptive mesh refinement/coarsening
              mesh_refinement->refine_and_coarsen_elements();
            }
          else // Refine based on reaching a target number of elements
            {
              //If we have reached the desired number of elements, we
              //dont do any more adaptive steps
              if (mesh.n_active_elem() >= param.nelem_target)
                {
                  libMesh::out << "We reached the target number of elements." << std::endl <<std::endl;
                  break;
                }

              mesh_refinement->flag_elements_by_nelem_target (error);

              // Carry out the adaptive mesh refinement/coarsening
              mesh_refinement->refine_and_coarsen_elements();
            }

          equation_systems.reinit();
        }
    }

  write_output_footers(param);

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;

}
