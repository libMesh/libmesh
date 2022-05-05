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

// <h1>Vector Finite Element Example 1 - Solving an uncoupled Poisson Problem using DG</h1>
// \author Alexander Lindsay
// \date 2019
//
// This is the fifth vector FE example program.  It solves the same PDE as the first vector example,
// but uses a vector monomial basis. The approximate solution of both ex1 and this example will
// converge to the same true solution, but at different rates. Moreover, although the PDE is linear,
// we demonstrate use of a non-linear solver in this example while ex1 uses an appropriate linear
// solver

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/parameters.h"
#include "libmesh/getpot.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/enum_elem_type.h"

namespace libMesh
{
template <typename>
class NumericVector;
template <typename>
class SparseMatrix;
}

// Bring in everything from the libMesh namespace
using namespace libMesh;

// residual assembly function
void compute_residual(const NumericVector<Number> & X,
                      NumericVector<Number> & R,
                      NonlinearImplicitSystem & S);

// jacobian assembly function
void compute_jacobian(const NumericVector<Number> & X,
                      SparseMatrix<Number> & J,
                      NonlinearImplicitSystem & S);

int
main(int argc, char ** argv)
{
  // Initialize libraries.
  LibMeshInit init(argc, argv);

  // This example requires a non-linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS ||
                               libMesh::default_solver_package() == TRILINOS_SOLVERS,
                           "--enable-petsc or --enable-trilinos");

  // Parse input file
  GetPot input_file("vector_fe_ex5.in");

  // But allow the command line to override it.
  input_file.parse_command_line(argc, argv);

  // Read DG parameters
  const Real epsilon = input_file("epsilon", -1);
  const Real sigma = input_file("sigma", 6);

  // Read mesh size
  const std::size_t nx = input_file("nx", 15);
  const std::size_t ny = input_file("ny", 15);

  // Read approximation order - default to FIRST
  const unsigned int approx_order = input_file("approx_order", 1);

  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i = 1; i < argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD4 elements. Note that at the end of this
  // fuction call, the mesh will be prepared and remote elements will be deleted
  MeshTools::Generation::build_square(mesh, nx, ny, -1., 1., -1., 1., QUAD4);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);
  equation_systems.parameters.set<Real>("epsilon") = epsilon;
  equation_systems.parameters.set<Real>("sigma") = sigma;

  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  auto & nl_system = equation_systems.add_system<NonlinearImplicitSystem>("Poisson");

  // Adds the variable "u" to "Poisson". Variable "u" will be approximated using
  // first-order (can be overridden with command line) vector monomial elements.
  // Since the mesh is 2-D, "u" will have two components.
  nl_system.add_variable("u", static_cast<Order>(approx_order), MONOMIAL_VEC);

  // Set the residual and Jacobian evaluation functions
  auto & nl_solver = *nl_system.nonlinear_solver;
  nl_solver.residual = compute_residual;
  nl_solver.jacobian = compute_jacobian;

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  nl_system.solve();

#ifdef LIBMESH_HAVE_EXODUS_API
  if (! approx_order) /* CONSTANT, MONOMIAL_VEC */
    {
      // Warn about trivial solution for CONSTANT approximations (Poisson must be at least C1)
      libMesh::out << "\nWARNING: The elemental vector order is CONSTANT. The solution" << std::endl
                   << "written out to 'out_constant.e' will be trivial." << std::endl;

      // Need to get string of variable names (vector components) for 'set_output_variables()'
      std::vector<std::string> varnames;
      equation_systems.build_variable_names(varnames);

      // Create ExodusII object with an unambigous filename (indicating this solution is special)
      ExodusII_IO exo_ptr(mesh);
      exo_ptr.write("out_constant.e");
      exo_ptr.set_output_variables(varnames);
      exo_ptr.write_element_data(equation_systems);
    }
  else
    ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
#endif

  // All done.
  return 0;
}
