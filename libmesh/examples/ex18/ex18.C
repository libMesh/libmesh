/* $Id: ex18.C,v 1.10 2006-07-20 22:48:46 roystgnr Exp $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

 // <h1>Example 18 - Unsteady Navier-Stokes Equations with DiffSystem</h1>
 //
 // This example shows how the transient nonlinear problem from
 // example 13 can be solved using the new (and experimental)
 // DiffSystem class framework

// Basic include files
#include "equation_systems.h"
#include "error_vector.h"
#include "getpot.h"
#include "gmv_io.h"
#include "kelly_error_estimator.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_refinement.h"
#include "uniform_refinement_estimator.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// The systems and solvers we may use
#include "naviersystem.h"
#include "diff_solver.h"
#include "euler_solver.h"
#include "steady_solver.h"

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  libMesh::init (argc, argv);
  {    
    // Parse the input file
    GetPot infile("ex18.in");

    // Read in parameters from the input file
    const Real global_tolerance          = infile("global_tolerance", 0.);
    const unsigned int nelem_target      = infile("n_elements", 400);
    const bool transient                 = infile("transient", true);
    const Real deltat                    = infile("deltat", 0.005);
    unsigned int n_timesteps             = infile("n_timesteps", 20);
    const unsigned int write_interval    = infile("write_interval", 5);
    const unsigned int coarsegridsize    = infile("coarsegridsize", 1);
    const unsigned int coarserefinements = infile("coarserefinements", 0);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 10);

    // Create a two-dimensional mesh.
    Mesh mesh (2);
    
    // And an object to refine it
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.coarsen_by_parents() = true;
    mesh_refinement.absolute_global_tolerance() = global_tolerance;
    mesh_refinement.nelem_target() = nelem_target;

    // Use the MeshTools::Generation mesh generator to create a uniform
    // grid on the square [-1,1]^D.  We instruct the mesh generator
    // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
    // elements in 3D.  Building these higher-order elements allows
    // us to use higher-order approximation, as in example 3.
    MeshTools::Generation::build_square (mesh,
                                         coarsegridsize,
                                         coarsegridsize,
                                         0., 1.,
                                         0., 1.,
                                         QUAD9);

    mesh_refinement.uniformly_refine(coarserefinements);

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system "Navier-Stokes" and its variables.
    NavierSystem & system = 
      equation_systems.add_system<NavierSystem> ("Navier-Stokes");

    // Solve this as a time-dependent or steady system
    if (transient)
      system.time_solver =
        AutoPtr<TimeSolver>(new EulerSolver(system));
    else
      {
        system.time_solver =
          AutoPtr<TimeSolver>(new SteadySolver(system));
        assert(n_timesteps == 1);
      }

    // Initialize the system
    equation_systems.init ();

    // Set the time stepping options
    system.deltat = deltat;

    // And the nonlinear solver options
    DiffSolver &solver = *system.time_solver->diff_solver;
    solver.quiet = infile("solver_quiet", true);
    solver.max_nonlinear_iterations =
      infile("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance =
      infile("relative_step_tolerance", 1.e-8);
    solver.relative_residual_tolerance =
      infile("relative_residual_tolerance", 1.e-9);

    // And the linear solver options
    solver.max_linear_iterations =
      infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
      infile("initial_linear_tolerance", 1.e-3);

    // Print information about the system to the screen.
    equation_systems.print_info();

    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
      {
        // A pretty update message
        std::cout << " Solving time step " << t_step << ", time = "
                  << system.time << std::endl;

        // Adaptively solve the timestep
        unsigned int a_step = 0;
        for (; a_step != max_adaptivesteps; ++a_step)
          {
            system.solve();

            ErrorVector error;

            AutoPtr<ErrorEstimator> error_estimator;

            // To solve to a tolerance in this problem we
            // need a better estimator than Kelly
            if (global_tolerance != 0.)
              {
                // We can't adapt to both a tolerance and a mesh
                // size at once
                assert (nelem_target == 0);

                UniformRefinementEstimator *u =
                  new UniformRefinementEstimator;

                // The lid-driven cavity problem isn't in H1, so
                // lets estimate H0 (i.e. L2) error
                u->sobolev_order() = 0;

                error_estimator.reset(u);
              }
            else
              {
                // If we aren't adapting to a tolerance we need a
                // target mesh size
                assert (nelem_target > 0);

                // Kelly is a lousy estimator to use for a problem
                // not in H1 - if we were doing more than a few
                // timesteps we'd need to turn off or limit the
                // maximum level of our adaptivity eventually
                error_estimator.reset(new KellyErrorEstimator);
              }

            // Calculate error based on u and v but not p
            error_estimator->component_scale.push_back(1.0); // u
            error_estimator->component_scale.push_back(1.0); // v
            error_estimator->component_scale.push_back(0.0); // p

            error_estimator->estimate_error(system, error);

            // Print out status at each adaptive step.
            Real global_error = error.l2_norm();
            std::cerr << "adaptive step " << a_step << ": ";
            if (global_tolerance != 0.)
              std::cerr << "global_error = " << global_error
                        << " with ";
            std::cerr << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
                      << " active dofs." << std::endl;
            if (global_tolerance != 0.)
              std::cerr << "worst element error = " << error.maximum()
                        << ", mean = " << error.mean() << std::endl;

            if (global_tolerance != 0.)
              {
                // If we've reached our desired tolerance, we
                // don't need any more adaptive steps
                if (global_error < global_tolerance)
                  break;
                mesh_refinement.flag_elements_by_error_tolerance(error);
              }
            else
              {
                // If flag_elements_by_nelem_target returns true, this
                // should be our last adaptive step.
                if (mesh_refinement.flag_elements_by_nelem_target(error))
                  {
                    mesh_refinement.refine_and_coarsen_elements();
                    equation_systems.reinit();
                    break;
                  }
              }

            // Carry out the adaptive mesh refinement/coarsening
            mesh_refinement.refine_and_coarsen_elements();
            equation_systems.reinit();
          }
        // Do one last solve if necessary
        if (a_step == max_adaptivesteps)
          {
            system.solve();
          }

        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();

        // Write out this timestep if we're requested to
        if ((t_step+1)%write_interval == 0)
          {
            OStringStream file_name;

            // We write the file name in the gmv auto-read format.
            file_name << "out.gmv.";
            OSSRealzeroright(file_name,3,0, t_step + 1);

            GMVIO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
          }
      }
  }
  
  // All done.  
  return libMesh::close ();
}
