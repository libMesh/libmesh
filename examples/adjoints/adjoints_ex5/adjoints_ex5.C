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



// Adjoints Example 5 - SolutionHistory, General Localized Vectors and Unsteady Adjoints.
// \author Vikram Garg
// \date 2013
//
// This example showcases three new capabilities in libMesh. The
// primary motivation is adjoint sensitivity analysis for unsteady
// problems. The PDE we are interested in is the simple 2-d heat
// equation:
// partial(T)/partial(t) - K Laplacian(T) = 0
// with initial condition:
// T(x,y;0) = sin(pi*x) sin(pi*y) and boundary conditions:
// T(boundary;t) = 0

// For these initial and boundary conditions, the exact solution
// u = exp(-K pi^2 t) * sin(pi*x) * sin(pi*y)

// We specify our Quantity of Interest (QoI) as
// Q(u) = int_{domain} u(x,y;1) sin(pi*x) sin(pi*y) dx dy, and
// are interested in computing the sensitivity dQ/dK

// For this QoI, the continuous adjoint problem reads,
// -partial(z)/partial(t) - K Laplacian(z) = 0
// with initial condition:
// T(x,y;1) = sin(pi*x) sin(pi*y)
// and boundary condition:
// T(boundary;t) = 0

// which has the exact solution,
// z = exp(-K pi^2 (1 - t)) * sin(pi*x) * sin(pi*y)
// which is the mirror image in time of the forward solution

// For an adjoint consistent space-time formulation, the discrete
// adjoint can be obtained by marching backwards from the adjoint
// initial condition and solving the transpose of the discrete primal
// problem at the last nonlinear solve of the corresponding primal
// timestep. This necessitates the storage of the primal solution at
// all timesteps, which is accomplished here using a
// MemorySolutionHistory object. As the name suggests, this object
// simply stores the primal solution (and other vectors we may choose
// to save), so that we can retrieve them later, whenever necessary.

// The discrete adjoint system for implicit time steppers requires the
// localization of vectors other than system.solution, which is
// accomplished using the localize_vectors method. In this particular
// example, we use the localized adjoint solution to assemble the
// residual contribution for the current adjoint timestep from the last
// computed adjoint timestep.

// Finally, The adjoint_advance_timestep method, the backwards time
// analog of advance_timestep prepares the time solver for solving the
// adjoint system, while the retrieve_timestep method retrieves the
// saved solutions at the current system.time, so that the adjoint
// sensitivity contribution for the current time can be computed.

// Local includes
#include "initial.h"
#include "adjoint_initial.h"
#include "femparameters.h"
#include "heatsystem.h"

// Libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/system_norm.h"
#include "libmesh/numeric_vector.h"

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

// SolutionHistory Includes
#include "libmesh/solution_history.h"
#include "libmesh/memory_solution_history.h"

// C++ includes
#include <iostream>
#include <sys/time.h>
#include <iomanip>

void write_output(EquationSystems & es,
                  unsigned int a_step,       // The adaptive step count
                  std::string solution_type) // primal or adjoint solve
{
#ifdef LIBMESH_HAVE_GMV
  MeshBase & mesh = es.get_mesh();

  std::ostringstream file_name_gmv;
  file_name_gmv << solution_type
                << ".out.gmv."
                << std::setw(2)
                << std::setfill('0')
                << std::right
                << a_step;

  GMVIO(mesh).write_equation_systems
    (file_name_gmv.str(), es);
#endif
}

void set_system_parameters(HeatSystem & system,
                           FEMParameters & param)
{
  // Use the prescribed FE type
  system.fe_family() = param.fe_family[0];
  system.fe_order() = param.fe_order[0];

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
      if (param.timesolver_core == "euler")
        {
          EulerSolver *eulersolver =
            new EulerSolver(system);

          eulersolver->theta = param.timesolver_theta;
          innersolver = eulersolver;
        }
      else
        libmesh_error_msg("This example (and unsteady adjoints in libMesh) only support Backward Euler and explicit methods.");

      system.time_solver =
        UniquePtr<TimeSolver>(innersolver);
    }
  else
    system.time_solver =
      UniquePtr<TimeSolver>(new SteadySolver(system));

  // The Memory Solution History object we will set the system SolutionHistory object to
  MemorySolutionHistory heatsystem_solution_history(system);
  system.time_solver->set_solution_history(heatsystem_solution_history);

  system.time_solver->reduce_deltat_on_diffsolver_failure =
    param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

  // Create any Dirichlet boundary conditions
  typedef
    std::map<boundary_id_type, FunctionBase<Number> *>::
    const_iterator Iter;

  for (Iter i = param.dirichlet_conditions.begin();
       i != param.dirichlet_conditions.end(); ++i)
    {
      boundary_id_type b = i->first;
      FunctionBase<Number> *f = i->second;
      std::set<boundary_id_type> bdys; bdys.insert(b);

      system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(bdys,
                                                                    param.dirichlet_condition_variables[b],
                                                                    f));

      libMesh::out << "Added Dirichlet boundary " << b << " for variables ";
      for (std::size_t vi=0; vi != param.dirichlet_condition_variables[b].size(); ++vi)
        libMesh::out << param.dirichlet_condition_variables[b][vi];
      libMesh::out << std::endl;
    }

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver() = UniquePtr<DiffSolver>(solver);
#else
      libmesh_error_msg("This example requires libMesh to be compiled with PETSc support.");
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver() = UniquePtr<DiffSolver>(solver);

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

// The main program.
int main (int argc, char** argv)
{
  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else
  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This doesn't converge with Trilinos for some reason...
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

  // Create a mesh with the given dimension, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm(), param.dimension);

  // And an object to refine it
  UniquePtr<MeshRefinement> mesh_refinement(new MeshRefinement(mesh));

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

  read_initial_parameters();

  system.project_solution(initial_value, initial_grad,
                          equation_systems.parameters);

  // Output the H1 norm of the initial conditions
  libMesh::out << "|U("
               << system.time
               << ")|= "
               << system.calculate_norm(*system.solution, 0, H1)
               << std::endl
               << std::endl;

  // Add an adjoint vector, this will be computed after the forward
  // time stepping is complete
  //
  // Tell the library not to save adjoint solutions during the forward
  // solve
  //
  // Tell the library not to project this vector, and hence, memory
  // solution history to not save it.
  //
  // Make this vector ghosted so we can localize it to each element
  // later.
  const std::string & adjoint_solution_name = "adjoint_solution0";
  system.add_vector("adjoint_solution0", false, GHOSTED);

  // Close up any resources initial.C needed
  finish_initialization();

  // Plot the initial conditions
  write_output(equation_systems, 0, "primal");

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

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
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          // Advance to the next timestep in a transient problem
          libMesh::out << "Advancing timestep" << std::endl << std::endl;
          system.time_solver->advance_timestep();

          // Write out this timestep
          write_output(equation_systems, t_step+1, "primal");
        }
      // End timestep loop

      ///////////////// Now for the Adjoint Solution //////////////////////////////////////

      // Now we will solve the backwards in time adjoint problem
      libMesh::out << std::endl << "Solving the adjoint problem" << std::endl;

      // We need to tell the library that it needs to project the adjoint, so
      // MemorySolutionHistory knows it has to save it

      // Tell the library to project the adjoint vector, and hence, memory solution history to
      // save it
      system.set_vector_preservation(adjoint_solution_name, true);

      libMesh::out << "Setting adjoint initial conditions Z("
                   << system.time
                   << ")"
                   <<std::endl;

      // Need to call adjoint_advance_timestep once for the initial condition setup
      libMesh::out<<"Retrieving solutions at time t="<<system.time<<std::endl;
      system.time_solver->adjoint_advance_timestep();

      // Output the H1 norm of the retrieved solutions (u^i and u^i+1)
      libMesh::out << "|U("
                   << system.time + system.deltat
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

      libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

      // The first thing we have to do is to apply the adjoint initial
      // condition. The user should supply these. Here they are specified
      // in the functions adjoint_initial_value and adjoint_initial_gradient
      system.project_vector(adjoint_initial_value,
                            adjoint_initial_grad,
                            equation_systems.parameters,
                            system.get_adjoint_solution(0));

      // Since we have specified an adjoint solution for the current
      // time (T), set the adjoint_already_solved boolean to true, so
      // we dont solve unneccesarily in the adjoint sensitivity method
      system.set_adjoint_already_solved(true);

      libMesh::out << "|Z("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(), 0, H1)
                   << std::endl
                   << std::endl;

      write_output(equation_systems, param.n_timesteps, "dual");

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

          // The adjoint_advance_timestep function calls the retrieve
          // function of the memory_solution_history class via the
          // memory_solution_history object we declared earlier.  The
          // retrieve function sets the system primal vectors to their
          // values at the current timestep.
          libMesh::out << "Retrieving solutions at time t=" << system.time << std::endl;
          system.time_solver->adjoint_advance_timestep();

          // Output the H1 norm of the retrieved solution
          libMesh::out << "|U("
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          libMesh::out << "|U("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                       << std::endl;

          system.set_adjoint_already_solved(false);

          system.adjoint_solve();

          // Now that we have solved the adjoint, set the
          // adjoint_already_solved boolean to true, so we dont solve
          // unneccesarily in the error estimator
          system.set_adjoint_already_solved(true);

          libMesh::out << "|Z("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(), 0, H1)
                       << std::endl
                       << std::endl;

          // Get a pointer to the primal solution vector
          NumericVector<Number> & primal_solution = *system.solution;

          // Get a pointer to the solution vector of the adjoint problem for QoI 0
          NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

          // Swap the primal and dual solutions so we can write out the adjoint solution
          primal_solution.swap(dual_solution_0);

          write_output(equation_systems, param.n_timesteps - (t_step + 1), "dual");

          // Swap back
          primal_solution.swap(dual_solution_0);
        }
      // End adjoint timestep loop

      // Now that we have computed both the primal and adjoint solutions, we compute the sensitivties to the parameter p
      // dQ/dp = partialQ/partialp - partialR/partialp
      // partialQ/partialp = (Q(p+dp) - Q(p-dp))/(2*dp), this is not supported by the library yet
      // partialR/partialp = (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp), where
      // R(u,z;p+dp) = int_{0}^{T} f(z;p+dp) - <partialu/partialt, z>(p+dp) - <g(u),z>(p+dp)
      // To do this we need to step forward in time, and compute the perturbed R at each time step and accumulate it
      // Then once all time steps are over, we can compute (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp)

      // Now we begin the timestep loop to compute the time-accurate
      // adjoint sensitivities
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          // A pretty update message
          libMesh::out << "Retrieving "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // Retrieve the primal and adjoint solutions at the current timestep
          system.time_solver->retrieve_timestep();

          libMesh::out << "|U("
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          libMesh::out << "|U("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                       << std::endl;

          libMesh::out << "|Z("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                       << std::endl
                       << std::endl;

          // Call the postprocess function which we have overloaded to compute
          // accumulate the perturbed residuals
          dynamic_cast<HeatSystem &>(system).perturb_accumulate_residuals(dynamic_cast<HeatSystem &>(system).get_parameter_vector());

          // Move the system time forward (retrieve_timestep does not do this)
          system.time += system.deltat;
        }

      // A pretty update message
      libMesh::out << "Retrieving final time = "
                   << system.time
                   << std::endl;

      // Retrieve the primal and adjoint solutions at the current timestep
      system.time_solver->retrieve_timestep();

      libMesh::out << "|U("
                   << system.time + system.deltat
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

      libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

      libMesh::out << "|Z("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                   << std::endl
                   << std::endl;

      // Call the postprocess function which we have overloaded to compute
      // accumulate the perturbed residuals
      dynamic_cast<HeatSystem &>(system).perturb_accumulate_residuals(dynamic_cast<HeatSystem &>(system).get_parameter_vector());

      // Now that we computed the accumulated, perturbed residuals, we can compute the
      // approximate sensitivity
      Number sensitivity_0_0 = (dynamic_cast<HeatSystem &>(system)).compute_final_sensitivity();

      // Print it out
      libMesh::out << "Sensitivity of QoI 0 w.r.t parameter 0 is: "
                   << sensitivity_0_0
                   << std::endl;

      // Hard coded assert to ensure that the actual numbers we are
      // getting are what they should be
      // The 2e-4 tolerance is chosen to ensure success even with
      // 32-bit floats
      libmesh_assert_less(std::abs(sensitivity_0_0 - (-5.37173)), 2.e-4);

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
