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



// <h1>FEMSystem Example 1 - Unsteady Navier-Stokes Equations with
// FEMSystem</h1>
// \author Roy Stogner
// \date 2006
//
// This example shows how the transient nonlinear problem from
// example 13 can be solved using the
// DifferentiableSystem class framework

// C++ includes
#include <iomanip>

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// The systems and solvers we may use
#include "naviersystem.h"
#include "libmesh/diff_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // This example fails without at least double precision FP
#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Trilinos and Eigen solvers NaN by default here.
  // We'll skip this example for now.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc");

  // Parse the input file
  GetPot infile("fem_system_ex1.in");

  // Read in parameters from the input file
  const Real global_tolerance          = infile("global_tolerance", 0.);
  const unsigned int nelem_target      = infile("n_elements", 400);
  const bool transient                 = infile("transient", true);
  const Real deltat                    = infile("deltat", 0.005);
  unsigned int n_timesteps             = infile("n_timesteps", 20);
  const unsigned int coarsegridsize    = infile("coarsegridsize", 1);
  const unsigned int coarserefinements = infile("coarserefinements", 0);
  const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 10);
  const unsigned int dim               = infile("dimension", 2);

#ifdef LIBMESH_HAVE_EXODUS_API
  const unsigned int write_interval    = infile("write_interval", 5);
#endif

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // We have only defined 2 and 3 dimensional problems
  libmesh_assert (dim == 2 || dim == 3);

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // And an object to refine it
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.coarsen_by_parents() = true;
  mesh_refinement.absolute_global_tolerance() = global_tolerance;
  mesh_refinement.nelem_target() = nelem_target;
  mesh_refinement.refine_fraction() = 0.3;
  mesh_refinement.coarsen_fraction() = 0.3;
  mesh_refinement.coarsen_threshold() = 0.1;

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements in 2D, or Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.
  if (dim == 2)
    MeshTools::Generation::build_square (mesh,
                                         coarsegridsize,
                                         coarsegridsize,
                                         0., 1.,
                                         0., 1.,
                                         QUAD9);
  else if (dim == 3)
    MeshTools::Generation::build_cube (mesh,
                                       coarsegridsize,
                                       coarsegridsize,
                                       coarsegridsize,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       HEX27);

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
    system.time_solver = libmesh_make_unique<EulerSolver>(system);
  else
    {
      system.time_solver = libmesh_make_unique<SteadySolver>(system);
      libmesh_assert_equal_to (n_timesteps, 1);
    }

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver & solver = *(system.time_solver->diff_solver().get());
  solver.quiet = infile("solver_quiet", true);
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations =
    infile("max_nonlinear_iterations", 15);
  solver.relative_step_tolerance =
    infile("relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance =
    infile("relative_residual_tolerance", 0.0);
  solver.absolute_residual_tolerance =
    infile("absolute_residual_tolerance", 0.0);

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
      libMesh::out << "\n\nSolving time step "
                   << t_step
                   << ", time = "
                   << system.time
                   << std::endl;

      // Adaptively solve the timestep
      unsigned int a_step = 0;
      for (; a_step != max_adaptivesteps; ++a_step)
        {
          system.solve();

          system.postprocess();

          ErrorVector error;

          std::unique_ptr<ErrorEstimator> error_estimator;

          // To solve to a tolerance in this problem we
          // need a better estimator than Kelly
          if (global_tolerance != 0.)
            {
              // We can't adapt to both a tolerance and a mesh
              // size at once
              libmesh_assert_equal_to (nelem_target, 0);

              UniformRefinementEstimator * u = new UniformRefinementEstimator;

              // The lid-driven cavity problem isn't in H1, so
              // lets estimate L2 error
              u->error_norm = L2;

              error_estimator.reset(u);
            }
          else
            {
              // If we aren't adapting to a tolerance we need a
              // target mesh size
              libmesh_assert_greater (nelem_target, 0);

              // Kelly is a lousy estimator to use for a problem
              // not in H1 - if we were doing more than a few
              // timesteps we'd need to turn off or limit the
              // maximum level of our adaptivity eventually
              error_estimator.reset(new KellyErrorEstimator);
            }

          // Calculate error based on u and v (and w?) but not p
          std::vector<Real> weights(2,1.0);  // u, v
          if (dim == 3)
            weights.push_back(1.0);          // w
          weights.push_back(0.0);            // p
          // Keep the same default norm type.
          std::vector<FEMNormType>
            norms(1, error_estimator->error_norm.type(0));
          error_estimator->error_norm = SystemNorm(norms, weights);

          error_estimator->estimate_error(system, error);

          // Print out status at each adaptive step.
          Real global_error = error.l2_norm();
          libMesh::out << "Adaptive step "
                       << a_step
                       << ": "
                       << std::endl;

          if (global_tolerance != 0.)
            libMesh::out << "Global_error = "
                         << global_error
                         << std::endl;

          if (global_tolerance != 0.)
            libMesh::out << "Worst element error = "
                         << error.maximum()
                         << ", mean = "
                         << error.mean()
                         << std::endl;

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
                  a_step = max_adaptivesteps;
                  break;
                }
            }

          // Carry out the adaptive mesh refinement/coarsening
          mesh_refinement.refine_and_coarsen_elements();
          equation_systems.reinit();

          libMesh::out << "Refined mesh to "
                       << mesh.n_active_elem()
                       << " active elements and "
                       << equation_systems.n_active_dofs()
                       << " active dofs."
                       << std::endl;
        }
      // Do one last solve if necessary
      if (a_step == max_adaptivesteps)
        {
          system.solve();

          system.postprocess();
        }

      // Advance to the next timestep in a transient problem
      system.time_solver->advance_timestep();

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out this timestep if we're requested to
      if ((t_step+1)%write_interval == 0)
        {
          std::ostringstream file_name;

          // We write the file in the ExodusII format.
          file_name << "out_"
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step+1
                    << ".e";

          ExodusII_IO(mesh).write_timestep(file_name.str(),
                                           equation_systems,
                                           1, // This number indicates how many time steps
                                              // are being written to the file
                                           system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
