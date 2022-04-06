// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Adjoints Example 7 - Unsteady Adjoint Refinement Error Estimator for Model Error Estimation.
// \author Vikram Garg
// \date 2020
//
// This example illustrates the use of the AdjointRefinementErrorEstimator in an unsteady setting. It also
// illustrates adjoint estimation of model error.
// The PDE we are interested in is the simple 2-d heat
// equation:
// partial(T)/partial(t) - K(x) Laplacian(T) = 1.0
// with initial condition:
// T(x,y;0) = 1.0 and boundary conditions:
// T(boundary;t) = 1.0
// with K(x) = 0.001 + x for the true model, and K = 0.01 for the coarse model.
// We are interested in estimating the effect of using K = 0.01 on 2 QoIs.

// We specify our Quantity of Interest (QoI) as
// Q0(u) = int_{domain} u(x,y;tf) 1.0 dx dy (A temporally non-smooth QoI evaluated at the final time)
// Q1(u) = int_{0}^{tf} int_{domain} u(x,y;t) 1.0 dx dy dt (A temporally smooth QoI)

// We verify the adjoint identity,
// Q(u_true) - Q(u_coarse) = R(u_true, z)

// Local includes
#include "initial.h"
#include "adjoint_initial.h"
#include "femparameters.h"
#include "heatsystem.h"
#include "sigma_physics.h"

// Libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/system_norm.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_solver_package.h"

#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/petsc_diff_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/twostep_time_solver.h"

#include "libmesh/getpot.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exodusII_io.h"

#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/error_vector.h"

// SolutionHistory Includes
#include "libmesh/solution_history.h"
#include "libmesh/memory_solution_history.h"
#include "libmesh/file_solution_history.h"

// C++ includes
#include <iostream>
#include <sys/time.h>
#include <iomanip>
#include <memory>

void write_output(EquationSystems & es,
                  unsigned int t_step,       // The current time step count
                  std::string solution_type, // primal or adjoint solve
                  FEMParameters & param)
{
  // Ignore parameters when there are no output formats available.
  libmesh_ignore(es, t_step, solution_type, param);

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      MeshBase & mesh = es.get_mesh();

      std::ostringstream file_name_gmv;
      file_name_gmv << solution_type
                    << ".out.gmv."
                    << std::setw(2)
                    << std::setfill('0')
                    << std::right
                    << t_step;

      GMVIO(mesh).write_equation_systems(file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_EXODUS_API
  if (param.output_exodus)
    {
      MeshBase & mesh = es.get_mesh();

      // We write out one file per timestep. The files are named in
      // the following way:
      // foo.e
      // foo.e-s002
      // foo.e-s003
      // ...
      // so that, if you open the first one with Paraview, it actually
      // opens the entire sequence of adapted files.
      std::ostringstream file_name_exodus;

      file_name_exodus << solution_type << ".e";
      if (t_step > 0)
        file_name_exodus << "-s"
                         << std::setw(3)
                         << std::setfill('0')
                         << std::right
                         << t_step + 1;

      // TODO: Get the current time from the System...
      ExodusII_IO(mesh).write_timestep(file_name_exodus.str(),
                                       es,
                                       1,
                                       /*time=*/t_step + 1);
    }
#endif
}

void set_system_parameters(HeatSystem &system, FEMParameters &param)
{
  // Use analytical jacobians?
  system.analytic_jacobians() = param.analytic_jacobians;

  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;
  system.numerical_jacobian_h = param.numerical_jacobian_h;

  // More desperate debugging options
  system.print_solution_norms    = param.print_solution_norms;
  system.print_solutions         = param.print_solutions;
  system.print_residual_norms    = param.print_residual_norms;
  system.print_residuals         = param.print_residuals;
  system.print_jacobian_norms    = param.print_jacobian_norms;
  system.print_jacobians         = param.print_jacobians;
  system.print_element_jacobians = param.print_element_jacobians;
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
        libmesh_error_msg("This example (and unsteady adjoints in libMesh) only support Backward Euler and explicit methods.");

      if (param.timesolver_tolerance)
        {
          TwostepTimeSolver *timesolver =
            new TwostepTimeSolver(system);

          timesolver->max_growth       = param.timesolver_maxgrowth;
          timesolver->target_tolerance = param.timesolver_tolerance;
          timesolver->upper_tolerance  = param.timesolver_upper_tolerance;
          timesolver->component_norm   = SystemNorm(param.timesolver_norm);

          timesolver->core_time_solver =
            std::unique_ptr<UnsteadySolver>(innersolver);

          system.time_solver =
            std::unique_ptr<UnsteadySolver>(timesolver);

          // The Memory/File Solution History object we will set the system SolutionHistory object to
          if(param.solution_history_type == "memory")
          {
            MemorySolutionHistory heatsystem_solution_history(system);
            system.time_solver->set_solution_history(heatsystem_solution_history);
          }
          else if (param.solution_history_type == "file")
          {
            FileSolutionHistory heatsystem_solution_history(system);
            system.time_solver->set_solution_history(heatsystem_solution_history);
          }
          else
          libmesh_error_msg("Unrecognized solution history type: " << param.solution_history_type);

        }
      else
      {
        system.time_solver =
          std::unique_ptr<TimeSolver>(innersolver);

        // The Memory/File Solution History object we will set the system SolutionHistory object to
        if(param.solution_history_type == "memory")
        {
          MemorySolutionHistory heatsystem_solution_history(system);
          system.time_solver->set_solution_history(heatsystem_solution_history);
        }
        else if (param.solution_history_type == "file")
        {
          FileSolutionHistory heatsystem_solution_history(system);
          system.time_solver->set_solution_history(heatsystem_solution_history);
        }
        else
        libmesh_error_msg("Unrecognized solution history type: " << param.solution_history_type);

      }
    }
  else
    system.time_solver =
      std::unique_ptr<TimeSolver>(new SteadySolver(system));

  system.time_solver->reduce_deltat_on_diffsolver_failure =
                                        param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

 #ifdef LIBMESH_ENABLE_DIRICHLET
  // Create any Dirichlet boundary conditions
  typedef
    std::map<boundary_id_type, FunctionBase<Number> *>::
    const_iterator Iter;

  for (Iter i = param.dirichlet_conditions.begin();
       i != param.dirichlet_conditions.end(); ++i)
    {
      boundary_id_type b = i->first;
      FunctionBase<Number> *f = i->second;

      system.get_dof_map().add_dirichlet_boundary(DirichletBoundary({b},
                                                                    param.dirichlet_condition_variables[b],
                                                                    f));

      libMesh::out << "Added Dirichlet boundary " << b << " for variables ";
      for (std::size_t vi=0; vi != param.dirichlet_condition_variables[b].size(); ++vi)
        libMesh::out << param.dirichlet_condition_variables[b][vi];
      libMesh::out << std::endl;
    }
 #endif // LIBMESH_ENABLE_DIRICHLET

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver() = std::unique_ptr<DiffSolver>(solver);
#else
      libmesh_error_msg("This example requires libMesh to be compiled with PETSc support.");
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver() = std::unique_ptr<DiffSolver>(solver);

      solver->quiet                       = param.solver_quiet;
      solver->verbose                     = param.solver_verbose;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
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
std::unique_ptr<AdjointRefinementEstimator>
build_adjoint_refinement_error_estimator(QoISet &qois, FEMPhysics* supplied_physics, FEMParameters &/*param*/)
{
  std::cout<<"Computing the error estimate using the Adjoint Refinement Error Estimator"<<std::endl<<std::endl;

  auto adjoint_refinement_estimator = std::make_unique<AdjointRefinementEstimator>();

  adjoint_refinement_estimator->qoi_set() = qois;

  // Set the residual evaluation physics pointer that belongs to the Adjoint Refinement
  // Error Estimator class to the residual_physics pointer passed by the user
  adjoint_refinement_estimator->set_residual_evaluation_physics(supplied_physics);

  // Also set the adjoint evaluation physics pointer to the same pointer
  adjoint_refinement_estimator->set_adjoint_evaluation_physics(supplied_physics);

  // We enrich the FE space for the dual problem by doing a uniform h refinements
  // Since we are dealing with model error, 0 h refinements are legal.
  adjoint_refinement_estimator->number_h_refinements = 0;

  return adjoint_refinement_estimator;
}
#endif // LIBMESH_ENABLE_AMR

// The main program.
int main (int argc, char ** argv)
{
  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_ignore(argc, argv);
  libmesh_example_requires(false, "--enable-amr");
#else
  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This doesn't converge with Trilinos for some reason...
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  libMesh::out << "Started " << argv[0] << std::endl;

  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("general.in");
    libmesh_error_msg_if(!i, '[' << init.comm().rank() << "] Can't find general.in; exiting early.");
  }
  GetPot infile("general.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the input file
  FEMParameters param(init.comm());
  param.read(infile);

  // Create a mesh with the given dimension, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm(), cast_int<unsigned char>(param.dimension));

  // And an object to refine it
  auto mesh_refinement = std::make_unique<MeshRefinement>(mesh);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  libMesh::out << "Building mesh" << std::endl;

  // Build a unit square
  ElemType elemtype;

  if (param.elementtype == "tri" ||
      param.elementtype == "unstructured")
    elemtype = TRI3;
  else
    elemtype = QUAD4;

  MeshTools::Generation::build_square (mesh, param.coarsegridx, param.coarsegridy,
                                       param.domain_xmin, param.domain_xmin + param.domain_edge_width,
                                       param.domain_ymin, param.domain_ymin + param.domain_edge_length,
                                       elemtype);

  libMesh::out << "Building system" << std::endl;

  HeatSystem & system = equation_systems.add_system<HeatSystem> ("HeatSystem");

  set_system_parameters(system, param);

  // Also build an FEMPhysics to provide the residual definition for the adjoint activities
  SigmaPhysics* sigma_physics = new SigmaPhysics();
  sigma_physics->init_data(system);

  libMesh::out << "Initializing systems" << std::endl;

  // Initialize the system
  equation_systems.init ();

  // Refine the grid again if requested
  for (unsigned int i=0; i != param.extrarefinements; ++i)
    {
      mesh_refinement->uniformly_refine(1);
      equation_systems.reinit();
    }

  libMesh::out << "Setting primal initial conditions" << std::endl;

  system.project_solution(initial_value, initial_grad,
                          equation_systems.parameters);

  // Output the H1 norm of the initial conditions
  libMesh::out << "|U("
               << system.time
               << ")|= "
               << system.calculate_norm(*system.solution, 0, H1)
               << std::endl
               << std::endl;

  // Add two adjoint vectors, these will be computed after the forward
  // and QoI calculation time stepping loops are complete
  //
  // Tell the library not to save adjoint solutions during the forward
  // solve
  //
  // Tell the library not to project this vector, and hence, memory/file
  // solution history to not save it during the forward run.
  //
  // Make this vector ghosted so we can localize it to each element
  // later.
  const std::string & adjoint_solution_name0 = "adjoint_solution0";
  const std::string & adjoint_solution_name1 = "adjoint_solution1";
  system.add_vector("adjoint_solution0", false, GHOSTED);
  system.add_vector("adjoint_solution1", false, GHOSTED);

  // To keep the number of vectors consistent between the primal and adjoint
  // loops, we will also pre-add the adjoint rhs vector
  system.add_vector("adjoint_rhs0", false, GHOSTED);
  system.add_vector("adjoint_rhs1", false, GHOSTED);

  // Plot the initial conditions
  write_output(equation_systems, 0, "primal", param);

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  // Accumulated integrated QoIs
  Number QoI_1_accumulated = 0.0;

  // In optimized mode we catch any solver errors, so that we can
  // write the proper footers before closing.  In debug mode we just
  // let the exception throw so that gdb can grab it.
#ifdef NDEBUG
  try
    {
#endif
      // Now we begin the timestep loop to compute the time-accurate
      // solution of the equations.
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          // A pretty update message
          libMesh::out << " Solving time step "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // Solve the forward problem at time t, to obtain the solution at time t + dt
          system.solve();

          // Output the H1 norm of the computed solution
          libMesh::out << "|U("
                      << system.time + system.time_solver->last_completed_timestep_size()
                      << ")|= "
                      << system.calculate_norm(*system.solution, 0, H1)
                      << std::endl;

          // Advance to the next timestep in a transient problem
          libMesh::out << "Advancing timestep" << std::endl << std::endl;
          system.time_solver->advance_timestep();

          // Write out this timestep
          write_output(equation_systems, t_step+1, "primal", param);
        }
      // End timestep loop

      // Set the time instant for the evaluation of the non-smooth QoI 0
      (dynamic_cast<HeatSystem &>(system).QoI_time_instant)[0] = system.time;

      // Now loop over the timesteps again to get the QoIs

      // Reset the initial time
      system.time = 0.0;

      system.time_solver->retrieve_timestep();

      QoI_1_accumulated = 0.0;

      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          system.time_solver->integrate_qoi_timestep();

          QoI_1_accumulated += (system.qoi)[1];
        }

      std::cout<< "The computed QoI 0 is " << std::setprecision(17) << (system.qoi)[0] << std::endl;
      std::cout<< "The computed QoI 1 is " << std::setprecision(17) << QoI_1_accumulated << std::endl;

      ///////////////// Now for the Adjoint Solution //////////////////////////////////////

      // Now we will solve the backwards in time adjoint problem
      libMesh::out << std::endl << "Solving the adjoint problems." << std::endl;

      // We need to tell the library that it needs to project the adjoint, so
      // MemorySolutionHistory knows it has to save it

      // Tell the library to project the adjoint vectors, and hence, solution history to
      // save them.
      system.set_vector_preservation(adjoint_solution_name0, true);
      system.set_vector_preservation(adjoint_solution_name1, true);

      libMesh::out << "Setting adjoint initial conditions Z("
                   << system.time
                   << ")"
                   <<std::endl;

      // The first thing we have to do is to apply the adjoint initial
      // condition. The user should supply these. Here they are specified
      // in the functions adjoint_initial_value and adjoint_initial_gradient
      system.project_vector(adjoint_initial_value0,
                            adjoint_initial_grad0,
                            equation_systems.parameters,
                            system.get_adjoint_solution(0));

      // The initial condition for adjoint 1 is zero, which is the default.

      // Since we have specified an adjoint solution for the current
      // time (T), set the adjoint_already_solved boolean to true, so
      // we dont solve unnecessarily in the adjoint sensitivity method
      system.set_adjoint_already_solved(true);

      libMesh::out << "|Z0("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                   << std::endl
                   << std::endl;

      libMesh::out << "|Z1("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(1), 0, H1)
                   << std::endl
                   << std::endl;

      // Call the adjoint advance timestep here to save the adjoint
      // initial conditions to disk
      system.time_solver->adjoint_advance_timestep();

      // Write out adjoint solutions for visualization.

      // Get a pointer to the primal solution vector
      NumericVector<Number> & primal_solution = *system.solution;

      // Get a pointer to the solution vector of the adjoint problem for QoI 0
      NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

      // Swap the primal and dual solutions so we can write out the adjoint solution
      primal_solution.swap(dual_solution_0);

      write_output(equation_systems, param.n_timesteps, "dual0", param);

      // Swap back
      primal_solution.swap(dual_solution_0);

      // Get a pointer to the solution vector of the adjoint problem for QoI 0
      NumericVector<Number> & dual_solution_1 = system.get_adjoint_solution(1);

      // Swap the primal and dual solutions so we can write out the adjoint solution
      primal_solution.swap(dual_solution_1);

      write_output(equation_systems, param.n_timesteps, "dual1", param);

      // Swap back
      primal_solution.swap(dual_solution_1);

      // Now that the adjoint initial condition is set, we will start the
      // backwards in time adjoint integration

      // For loop stepping backwards in time
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          //A pretty update message
          libMesh::out << " Solving adjoint time step "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // Output the H1 norm of the retrieved primal solution from the last call
          // to adjoint_advance_timestep
          libMesh::out << "|U("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          libMesh::out << "|U("
                         << system.time - (system.time_solver->last_completed_timestep_size())/((param.timesolver_tolerance) ? 2.0 : 1.0)
                         << ")|= "
                         << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                         << std::endl;

          system.set_adjoint_already_solved(false);

          // Swap in the physics we want to be used for the adjoint evaluation
          DifferentiablePhysics * adjoint_evaluation_physics = dynamic_cast<DifferentiablePhysics *>(sigma_physics);
          dynamic_cast<DifferentiableSystem &>(system).swap_physics(adjoint_evaluation_physics);

          system.adjoint_solve();

          // Swap back the physics
          dynamic_cast<DifferentiableSystem &>(system).swap_physics(adjoint_evaluation_physics);

          // Now that we have solved the adjoint, set the
          // adjoint_already_solved boolean to true, so we dont solve
          // unnecessarily in the error estimator
          system.set_adjoint_already_solved(true);

          libMesh::out << "Saving adjoint and retrieving primal solutions at time t=" << system.time - system.deltat << std::endl;

          // The adjoint_advance_timestep function calls the retrieve and store
          // function of the memory_solution_history class via the
          // memory_solution_history object we declared earlier.  The
          // retrieve function sets the system primal vectors to their
          // values at the current time instance.
          system.time_solver->adjoint_advance_timestep();

          libMesh::out << "|Z0("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                       << std::endl
                       << std::endl;

          libMesh::out << "|Z1("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(1), 0, H1)
                       << std::endl
                       << std::endl;

          // Get a pointer to the primal solution vector
          primal_solution = *system.solution;

          // Get a pointer to the solution vector of the adjoint problem for QoI 0
          dual_solution_0 = system.get_adjoint_solution(0);

          // Swap the primal and dual solutions so we can write out the adjoint solution
          primal_solution.swap(dual_solution_0);

          write_output(equation_systems, param.n_timesteps - (t_step + 1), "dual0", param);

          // Swap back
          primal_solution.swap(dual_solution_0);

          // Get a pointer to the solution vector of the adjoint problem for QoI 1
          dual_solution_1 = system.get_adjoint_solution(1);

          // Swap the primal and dual solutions so we can write out the adjoint solution
          primal_solution.swap(dual_solution_1);

          write_output(equation_systems, param.n_timesteps - (t_step + 1), "dual1", param);

          // Swap back
          primal_solution.swap(dual_solution_1);
        }
      // End adjoint timestep loop

      // Reset the time before looping over time again
      system.time = 0.0;

      // Now that we have computed both the primal and adjoint solutions, we can compute the goal-oriented error estimates.
      // For this, we will need to build a ARefEE error estimator object, and supply a pointer to the 'true physics' object,
      // as we did for the adjoint solve. This object will now define the residual for the dual weighted error evaluation at
      // each timestep.
      // For adjoint consistency, it is important that we maintain the same integration method for the error estimation integral
      // as we did for the primal and QoI integration. This is ensured within the library, but is currently ONLY supported for
      // the Backward-Euler (theta = 0.5) time integration method.
      QoISet qois;

      qois.add_indices({0,1});
      qois.set_weight(0, 0.5);
      qois.set_weight(1, 0.5);

      // The adjoint refinement error estimator object (which also computed model error)
      auto adjoint_refinement_error_estimator =
        build_adjoint_refinement_error_estimator(qois, dynamic_cast<FEMPhysics *>(sigma_physics), param);

      // Error Vector to be filled by the error estimator at each timestep
      ErrorVector QoI_elementwise_error;
      std::vector<Number> QoI_spatially_integrated_error(system.n_qois());

      // Error Vector to hold the total accumulated error
      ErrorVector accumulated_QoI_elementwise_error;
      std::vector<Number> accumulated_QoI_spatially_integrated_error(system.n_qois());

      // Retrieve the primal and adjoint solutions at the current timestep
      system.time_solver->retrieve_timestep();

      // A pretty update message
      libMesh::out << "Retrieved, "
               << "time = "
               << system.time
               << std::endl;

      libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

      libMesh::out << "|U_old("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

      libMesh::out << "|Z0("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                   << std::endl
                   << std::endl;

      libMesh::out << "|Z1("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(1), 0, H1)
                   << std::endl
                   << std::endl;

      // Now we begin the timestep loop to compute the time-accurate
      // adjoint sensitivities
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          system.time_solver->integrate_adjoint_refinement_error_estimate(*adjoint_refinement_error_estimator, QoI_elementwise_error);

          // A pretty update message
          libMesh::out << "Retrieved, "
               << "time = "
               << system.time
               << std::endl;

          libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

          libMesh::out << "|U_old("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

          libMesh::out << "|Z0("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                   << std::endl
                   << std::endl;

          libMesh::out << "|Z1("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(1), 0, H1)
                   << std::endl
                   << std::endl;

          // Error contribution from this timestep
          for(unsigned int i = 0; i < QoI_elementwise_error.size(); i++)
            accumulated_QoI_elementwise_error[i] += QoI_elementwise_error[i];

          // QoI wise error contribution from this timestep
          for (auto j : make_range(system.n_qois()))
          {
            // Skip this QoI if not in the QoI Set
            if ((adjoint_refinement_error_estimator->qoi_set()).has_index(j))
            {
              accumulated_QoI_spatially_integrated_error[j] += (system.qoi_error_estimates)[j];
            }
          }
        }

      std::cout<<"Time integrated error estimate for QoI 0: "<<std::setprecision(17)<<accumulated_QoI_spatially_integrated_error[0]<<std::endl;
      std::cout<<"Time integrated error estimate for QoI 1: "<<std::setprecision(17)<<accumulated_QoI_spatially_integrated_error[1]<<std::endl;

      // Hard coded test to ensure that the actual numbers we are
      // getting are what they should be
      // The 2e-4 tolerance is chosen to ensure success even with
      // 32-bit floats
      if(param.timesolver_tolerance)
      {
        libmesh_error_msg_if(std::abs(system.time - (1.7548735069535084)) >= 2.e-4,
                             "Mismatch in end time reached by adaptive timestepper!");

        libmesh_error_msg_if(std::abs(accumulated_QoI_spatially_integrated_error[0] - (-0.75139371165754232)) >= 2.e-4,
                             "Error Estimator identity not satisfied!");

        libmesh_error_msg_if(std::abs(accumulated_QoI_spatially_integrated_error[1] - (-0.70905030091469889)) >= 2.e-4,
                             "Error Estimator identity not satisfied!");
      }
      else
      {
        libmesh_error_msg_if(std::abs(accumulated_QoI_spatially_integrated_error[0] - (-0.44809323880671958)) >= 2.e-4,
                             "Error Estimator identity not satisfied!");

        libmesh_error_msg_if(std::abs(accumulated_QoI_spatially_integrated_error[1] - (-0.23895256319278346)) >= 2.e-4,
                             "Error Estimator identity not satisfied!");
      }
#ifdef NDEBUG
    }
  catch (...)
    {
      libMesh::err << '[' << mesh.processor_id()
                   << "] Caught exception; exiting early." << std::endl;
    }
#endif

  libMesh::err << '[' << mesh.processor_id()
               << "] Completing output."
               << std::endl;

  // All done.
  return 0;

#endif // LIBMESH_ENABLE_AMR
}
