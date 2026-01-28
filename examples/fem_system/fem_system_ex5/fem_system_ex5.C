// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>FEMSystem Example 5 - Unsteady Linear Elasticity with
// IsoGeometric Analysis</h1>
// \author Roy Stogner
// \date 2019
//
// This example shows how to solve the three-dimensional transient
// linear elasticity equations with elements suitable for IsoGeometric
// Analysis, as read from a DynaIO MeshInput.
//
// This is just FEMSystem Example 3 on an IGA mesh.

// C++ includes
#include <iomanip>
#include <memory>

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/dyna_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/numeric_vector.h"

// The systems and solvers we may use
#include "elasticity_system.h"
#include "libmesh/diff_solver.h"
#include "libmesh/newmark_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/elem.h"
#include "libmesh/newton_solver.h"
#include "libmesh/eigen_sparse_linear_solver.h"

#define x_scaling 1.3

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

  // This example requires 3D calculations
  libmesh_example_requires(LIBMESH_DIM > 2, "3D support");

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Parse the input file
  GetPot infile("fem_system_ex5.in");

  // Override input file arguments from the command line
  infile.parse_command_line(argc, argv);

  // Read in parameters from the input file
  const Real deltat           = infile("deltat", 0.25);
  unsigned int n_timesteps    = infile("n_timesteps", 1);

#ifdef LIBMESH_HAVE_EXODUS_API
  const unsigned int write_interval    = infile("write_interval", 1);
#endif

  // Initialize the mesh
  const unsigned int dim = infile("dim", 2);
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // Create a mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);

  // Use an IsoGeometricAnalysis mesh?
  const bool use_iga = infile("use_iga", true);

  // Declare a default FEType here to override if needed for IGA later
  FEType fe_type;

  const std::string iga_filename = infile("iga_filename","DIE!");

  // Load an IGA pressurized cylinder mesh or build a cantilever mesh
  if (use_iga)
    {
      DynaIO dyna_io(mesh);

      libMesh::out << "\nReading IGA mesh.\n" << std::endl;
      if (mesh.processor_id() == 0)
        dyna_io.read(iga_filename);
      MeshCommunication().broadcast (mesh);
      mesh.prepare_for_use();

      fe_type.order = 2;
      fe_type.family = RATIONAL_BERNSTEIN;
    }
  else
    {
      if (dim == 3)
        MeshTools::Generation::build_cube
          (mesh, 32, 8, 4, 0., 1.*x_scaling, 0., 0.3, 0., 0.1, HEX8);
      else if (dim == 2)
        MeshTools::Generation::build_square
          (mesh, 8, 4, 0., 1.*x_scaling, 0., 0.3, QUAD4);
      else
        libmesh_error_msg("Unsupported dim = " << dim);
    }

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Change our boundary id definitions to match the problem mesh
  if (dim == 2)
    {
      boundary_id_min_z = 99;
      boundary_id_min_y = 0;
      boundary_id_max_x = 1;
      boundary_id_max_y = 2;
      boundary_id_min_x = 3;
      boundary_id_max_z = 99;
    }

  // Let's add some node and edge boundary conditions.
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  for (const auto & elem : mesh.element_ptr_range())
    {
      if (!use_iga)
        {
          unsigned int side_max_x = 0, side_max_y = 0, side_min_z = 0;
          bool
            found_side_max_x = false, found_side_max_y = false,
            found_side_min_z = false;
          for (auto side : elem->side_index_range())
            {
              if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_x))
                {
                  side_max_x = side;
                  found_side_max_x = true;
                }

              if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_y))
                {
                  side_max_y = side;
                  found_side_max_y = true;
                }

              if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_min_z)
                  || dim == 2)
                {
                  side_min_z = side;
                  found_side_min_z = true;
                }
            }

          // If elem has sides on boundaries
          // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MIN_Z
          // then let's set a node boundary condition
          if (found_side_max_x && found_side_max_y && found_side_min_z)
            for (auto n : elem->node_index_range())
              if (elem->is_node_on_side(n, side_max_x) &&
                  elem->is_node_on_side(n, side_max_y) &&
                  (dim == 2 || elem->is_node_on_side(n, side_min_z)))
                mesh.get_boundary_info().add_node(elem->node_ptr(n), node_boundary_id);

          // If elem has sides on boundaries
          // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Z
          // then let's set an edge boundary condition
          if (found_side_max_x && found_side_min_z)
            for (auto e : elem->edge_index_range())
              if (elem->is_edge_on_side(e, side_max_x) &&
                  (dim == 2 || elem->is_edge_on_side(e, side_min_z)))
                mesh.get_boundary_info().add_edge(elem, e, edge_boundary_id);
        }
      else // if (use_iga)
        {
          // On our IGA mesh, we don't have BCs set yet, but:
          // Side 0 is the r=1 outer quarter circle
          // Side 1 is the x=0 symmetry boundary
          // Side 2 is the r=0.5 inner quarter circle
          // Side 3 is the y=0 symmetry boundary

          // The bottom side of our IGA mesh file should have v
          // displacement fixed
          if (elem->type() == QUAD9 && !elem->neighbor_ptr(3))
            mesh.get_boundary_info().add_side(elem, 3, fixed_v_boundary_id);

          // The left/top side of our IGA mesh file should have u
          // displacement fixed
          if (elem->type() == QUAD9 && !elem->neighbor_ptr(1))
            mesh.get_boundary_info().add_side(elem, 1, fixed_u_boundary_id);

          // The inside of our cylinder should have pressure applied
          if (elem->type() == QUAD9 && !elem->neighbor_ptr(2))
            mesh.get_boundary_info().add_side(elem, 2, pressure_boundary_id);
        }
    }

  // We're done modifying the BoundaryInfo; get its caches up to date.
  mesh.get_boundary_info().regenerate_id_sets();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system "Linear Elasticity" and its variables.
  ElasticitySystem & system =
    equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");

  system.set_dim(dim);

  system.set_fe_type(fe_type);

  // Solve this as a time-dependent or steady system
  std::string time_solver = infile("time_solver","DIE!");

  ExplicitSystem * v_system = nullptr;
  ExplicitSystem * a_system = nullptr;

  if (time_solver == "newmark")
    {
      // Create ExplicitSystem to help output velocity
      v_system = &equation_systems.add_system<ExplicitSystem> ("Velocity");
      v_system->add_variable("u_vel", fe_type);
      v_system->add_variable("v_vel", fe_type);
      if (dim == 3)
        v_system->add_variable("w_vel", fe_type);

      // Create ExplicitSystem to help output acceleration
      a_system = &equation_systems.add_system<ExplicitSystem> ("Acceleration");
      a_system->add_variable("u_accel", fe_type);
      a_system->add_variable("v_accel", fe_type);
      if (dim == 3)
        a_system->add_variable("w_accel", fe_type);
    }

  if (time_solver == "newmark")
    system.time_solver = std::make_unique<NewmarkSolver>(system);

  else if (time_solver == "euler")
    {
      system.time_solver = std::make_unique<EulerSolver>(system);
      EulerSolver & euler_solver = cast_ref<EulerSolver &>(*(system.time_solver.get()));
      euler_solver.theta = infile("theta", 1.0);
    }

  else if (time_solver == "euler2")
    {
      system.time_solver = std::make_unique<Euler2Solver>(system);
      Euler2Solver & euler_solver = cast_ref<Euler2Solver &>(*(system.time_solver.get()));
      euler_solver.theta = infile("theta", 1.0);
    }

  else if (time_solver == "steady")
    {
      system.time_solver = std::make_unique<SteadySolver>(system);
      libmesh_assert_equal_to (n_timesteps, 1);
    }
  else
    libmesh_error_msg(std::string("ERROR: invalid time_solver ")+time_solver);

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver & solver = *(system.time_solver->diff_solver().get());
  solver.quiet = infile("solver_quiet", true);
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations = infile("max_nonlinear_iterations", 15);
  solver.relative_step_tolerance = infile("relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance = infile("relative_residual_tolerance", 0.0);
  solver.absolute_residual_tolerance = infile("absolute_residual_tolerance", 0.0);

  // And the linear solver options
  solver.max_linear_iterations = infile("max_linear_iterations", 50000);
  solver.initial_linear_tolerance = infile("initial_linear_tolerance", 1.e-3);

  // And the system
  system.print_element_jacobians = infile("print_element_jacobians", false);
  system.print_element_solutions = infile("print_element_solutions", false);
  system.print_element_residuals = infile("print_element_residuals", false);
  system.print_residuals = infile("print_residuals", false);
  system.print_jacobians = infile("print_residuals", false);
  system.print_residual_norms = infile("print_residual_norms", false);
  system.print_jacobian_norms = infile("print_residuals", false);

  // Print information about the system to the screen.
  equation_systems.print_info();

  // If we're using EulerSolver or Euler2Solver, and we're using EigenSparseLinearSolver,
  // then we need to reset the EigenSparseLinearSolver to use SPARSELU because BICGSTAB
  // chokes on the Jacobian matrix we give it and Eigen's GMRES currently doesn't work.
  NewtonSolver * newton_solver = dynamic_cast<NewtonSolver *>( &solver );
  if (newton_solver &&
      (time_solver == "euler" || time_solver == "euler2"))
    {
#ifdef LIBMESH_HAVE_EIGEN_SPARSE
      LinearSolver<Number> & linear_solver = newton_solver->get_linear_solver();
      EigenSparseLinearSolver<Number> * eigen_linear_solver =
        dynamic_cast<EigenSparseLinearSolver<Number> *>(&linear_solver);

      if (eigen_linear_solver )
        eigen_linear_solver->set_solver_type(SPARSELU);
#endif
    }

  if (time_solver == "newmark")
    {
      NewmarkSolver * newmark = cast_ptr<NewmarkSolver*>(system.time_solver.get());
      newmark->compute_initial_accel();

      // Copy over initial velocity and acceleration for output.
      // Note we can do this because of the matching variables/FE spaces
      *(v_system->solution) = system.get_vector("_old_solution_rate");
      *(a_system->solution) = system.get_vector("_old_solution_accel");
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  // Output initial state (if we can - IGA meshs don't serialize yet)
  if (mesh.is_serial() || !use_iga)
    {
      std::ostringstream file_name;

      // We write the file in the ExodusII format.
      file_name << std::string("out.")+time_solver+std::string(".e-s.")
                << std::setw(3)
                << std::setfill('0')
                << std::right
                << 0;

      ExodusII_IO(mesh).write_timestep(file_name.str(),
                                       equation_systems,
                                       1, // This number indicates how many time steps
                                          // are being written to the file
                                       system.time);
    }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

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

      system.solve();

      // Advance to the next timestep in a transient problem
      system.time_solver->advance_timestep();

      // Copy over updated velocity and acceleration for output.
      // Note we can do this because of the matching variables/FE spaces
      if (time_solver == "newmark")
        {
          *(v_system->solution) = system.get_vector("_old_solution_rate");
          *(a_system->solution) = system.get_vector("_old_solution_accel");
        }

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out this timestep if we're requested to
      // (if we can - IGA meshs don't serialize yet)
      if ((mesh.is_serial() || !use_iga) &&
          ((t_step+1)%write_interval == 0))
        {
          std::ostringstream file_name;

          // We write the file in the ExodusII format.
          file_name << std::string("out.")+time_solver+std::string(".e-s.")
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step+1;

          ExodusII_IO(mesh).write_timestep(file_name.str(),
                                           equation_systems,
                                           1, // This number indicates how many time steps
                                              // are being written to the file
                                           system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }

  // All done.
  return 0;
}
