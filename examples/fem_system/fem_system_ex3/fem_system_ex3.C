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



// <h1>FEMSystem Example 3 - Unsteady Linear Elasticity with
// FEMSystem</h1>
// \author Paul Bauman
// \date 2015
//
// This example shows how to solve the three-dimensional transient
// linear elasticity equations using the DifferentiableSystem class framework.
// This is just Systems of Equations Example 6 recast.

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

  // Parse the input file
  GetPot infile("fem_system_ex3.in");

  // Override input file arguments from the commannd line
  infile.parse_command_line(argc, argv);

  // Read in parameters from the input file
  const Real deltat           = infile("deltat", 0.25);
  unsigned int n_timesteps    = infile("n_timesteps", 1);

#ifdef LIBMESH_HAVE_EXODUS_API
  const unsigned int write_interval    = infile("write_interval", 1);
#endif

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     32,
                                     8,
                                     4,
                                     0., 1.*x_scaling,
                                     0., 0.3,
                                     0., 0.1,
                                     HEX8);


  // Print information about the mesh to the screen.
  mesh.print_info();

  // Let's add some node and edge boundary conditions.
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      unsigned int
        side_max_x = 0, side_min_y = 0,
        side_max_y = 0, side_max_z = 0;
      bool
        found_side_max_x = false, found_side_max_y = false,
        found_side_min_y = false, found_side_max_z = false;
      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
            {
              side_max_x = side;
              found_side_max_x = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MIN_Y))
            {
              side_min_y = side;
              found_side_min_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Y))
            {
              side_max_y = side;
              found_side_max_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }
        }

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z
      // then let's set a node boundary condition
      if (found_side_max_x && found_side_max_y && found_side_max_z)
        for (unsigned int n=0; n<elem->n_nodes(); n++)
          if (elem->is_node_on_side(n, side_max_x) &&
              elem->is_node_on_side(n, side_max_y) &&
              elem->is_node_on_side(n, side_max_z))
            mesh.get_boundary_info().add_node(elem->node_ptr(n), NODE_BOUNDARY_ID);

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
      // then let's set an edge boundary condition
      if (found_side_max_x && found_side_min_y)
        for (unsigned int e=0; e<elem->n_edges(); e++)
          if (elem->is_edge_on_side(e, side_max_x) &&
              elem->is_edge_on_side(e, side_min_y))
            mesh.get_boundary_info().add_edge(elem, e, EDGE_BOUNDARY_ID);
    }

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system "Navier-Stokes" and its variables.
  ElasticitySystem & system =
    equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");

  // Solve this as a time-dependent or steady system
  std::string time_solver = infile("time_solver","DIE!");

  ExplicitSystem * v_system;
  ExplicitSystem * a_system;

  if( time_solver == std::string("newmark") )
    {
      // Create ExplicitSystem to help output velocity
      v_system = &equation_systems.add_system<ExplicitSystem> ("Velocity");
      v_system->add_variable("u_vel", FIRST, LAGRANGE);
      v_system->add_variable("v_vel", FIRST, LAGRANGE);
      v_system->add_variable("w_vel", FIRST, LAGRANGE);

      // Create ExplicitSystem to help output acceleration
      a_system = &equation_systems.add_system<ExplicitSystem> ("Acceleration");
      a_system->add_variable("u_accel", FIRST, LAGRANGE);
      a_system->add_variable("v_accel", FIRST, LAGRANGE);
      a_system->add_variable("w_accel", FIRST, LAGRANGE);
    }

  if( time_solver == std::string("newmark"))
    system.time_solver.reset(new NewmarkSolver(system));

  else if( time_solver == std::string("euler") )
    {
      system.time_solver.reset(new EulerSolver(system));
      EulerSolver & euler_solver = libmesh_cast_ref<EulerSolver &>(*(system.time_solver.get()));
      euler_solver.theta = infile("theta", 1.0);
    }

  else if( time_solver == std::string("euler2") )
    {
      system.time_solver.reset(new Euler2Solver(system));
      Euler2Solver & euler_solver = libmesh_cast_ref<Euler2Solver &>(*(system.time_solver.get()));
      euler_solver.theta = infile("theta", 1.0);
    }

  else if( time_solver == std::string("steady"))
    {
      system.time_solver.reset(new SteadySolver(system));
      libmesh_assert_equal_to (n_timesteps, 1);
    }
  else
    libmesh_error_msg(std::string("ERROR: invalud time_solver ")+time_solver);

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

  // Print information about the system to the screen.
  equation_systems.print_info();

  // If we're using EulerSolver or Euler2Solver, and we're using EigenSparseLinearSolver,
  // then we need to reset the EigenSparseLinearSolver to use SPARSELU because BICGSTAB
  // chokes on the Jacobian matrix we give it and Eigen's GMRES currently doesn't work.
  NewtonSolver * newton_solver = dynamic_cast<NewtonSolver *>( &solver );
  if( newton_solver &&
      (time_solver == std::string("euler") || time_solver == std::string("euler2") ) )
    {
      LinearSolver<Number> & linear_solver = newton_solver->get_linear_solver();
      EigenSparseLinearSolver<Number> * eigen_linear_solver =
        dynamic_cast<EigenSparseLinearSolver<Number> *>(&linear_solver);

      if( eigen_linear_solver )
        eigen_linear_solver->set_solver_type(SPARSELU);
    }

  if( time_solver == std::string("newmark") )
    {
      NewmarkSolver * newmark = cast_ptr<NewmarkSolver*>(system.time_solver.get());
      newmark->compute_initial_accel();

      // Copy over initial velocity and acceleration for output.
      // Note we can do this because of the matching variables/FE spaces
      *(v_system->solution) = system.get_vector("_old_solution_rate");
      *(a_system->solution) = system.get_vector("_old_solution_accel");
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  // Output initial state
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
      if( time_solver == std::string("newmark") )
        {
          *(v_system->solution) = system.get_vector("_old_solution_rate");
          *(a_system->solution) = system.get_vector("_old_solution_accel");
        }

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out this timestep if we're requested to
      if ((t_step+1)%write_interval == 0)
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
