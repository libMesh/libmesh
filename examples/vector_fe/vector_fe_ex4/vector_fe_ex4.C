/* The libMesh Finite Element Library. */
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

 // <h1>Vector Finite Elements Example 3 - Nedelec Elements</h1>
 //
 // This example shows an example of using the Nedelec elements of the
 // first type to solve a model problem in H(curl).

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exact_solution.h"
#include "libmesh/string_to_enum.h"

// The systems and solvers we may use
#include "curl_curl_system.h"
#include "libmesh/diff_solver.h"
#include "libmesh/steady_solver.h"

#include "solution_function.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Parse the input file
  GetPot infile("vector_fe_ex4.in");

  // Read in parameters from the input file
  const unsigned int grid_size = infile( "grid_size", 2 );

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_assert(3 <= LIBMESH_DIM, "2D/3D support");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D. We must use TRI6 elements for the
  // Nedelec triangle elements.

  std::string elem_str = command_line_value( std::string("element_type"),
					     std::string("HEX27") );

  if( elem_str != "HEX20" && elem_str != "HEX27" )
    {
      libMesh::err << "This example must be run with HEX20 or HEX27."
		   << std::endl
		   << "You entered: " << elem_str << std::endl;
      libmesh_error();
    }

  MeshTools::Generation::build_cube (mesh,
                                     grid_size,
                                     grid_size,
                                     grid_size,
                                     -1., 1,
                                     -1., 1.,
                                     -1., 1,
                                     Utility::string_to_enum<ElemType>(elem_str));


  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system "Navier-Stokes" and its variables.
  CurlCurlSystem & system =
    equation_systems.add_system<CurlCurlSystem> ("CurlCurl");

  // This example only implements the steady-state problem
  system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(system));

  // Initialize the system
  equation_systems.init();

  // And the nonlinear solver options
  DiffSolver &solver = *(system.time_solver->diff_solver().get());
  solver.quiet = infile("solver_quiet", true);
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations =
    infile("max_nonlinear_iterations", 15);
  solver.relative_step_tolerance =
    infile("relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance =
    infile("relative_residual_tolerance", 1.0e-13);
  solver.absolute_residual_tolerance =
    infile("absolute_residual_tolerance", 0.0);

  // And the linear solver options
  solver.max_linear_iterations =
    infile("max_linear_iterations", 50000);
  solver.initial_linear_tolerance =
    infile("initial_linear_tolerance", 1.e-10);

  // Print information about the system to the screen.
  equation_systems.print_info();

  system.solve();

  ExactSolution exact_sol( equation_systems );

  std::vector<FunctionBase<Number>* > sols;
  std::vector<FunctionBase<Gradient>* > grads;

  sols.push_back( new SolutionFunction(system.variable_number("u")) );
  grads.push_back( new SolutionGradient(system.variable_number("u")) );

  exact_sol.attach_exact_values(sols);
  exact_sol.attach_exact_derivs(grads);

  // Use higher quadrature order for more accurate error results
  int extra_error_quadrature = infile("extra_error_quadrature",2);
  exact_sol.extra_quadrature_order(extra_error_quadrature);

  // Compute the error.
  exact_sol.compute_error("CurlCurl", "u");

  // Print out the error values
  std::cout << "L2-Error is: "
	    << exact_sol.l2_error("CurlCurl", "u")
	    << std::endl;
  std::cout << "HCurl semi-norm error is: "
	    << exact_sol.error_norm("CurlCurl", "u", HCURL_SEMINORM )
	    << std::endl;
  std::cout << "HCurl-Error is: "
	    << exact_sol.hcurl_error("CurlCurl", "u")
	    << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}
