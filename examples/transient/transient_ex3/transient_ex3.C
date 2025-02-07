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

// <h1>Transient Example 3 - DG/FV formulation of the 2D Advection Equation with Explicit timestepping</h1>
// \author David Knezevic
// \date 2012
// \author John W. Peterson
// \date 2025 (modernization and libmesh example)
//
// This example program demonstrates one way to implement an explicit
// timestepping Discontinous Galerkin/Finite Volume formulation of the
// 2D advection equation. The example comes with an input file which
// is set up as a finite volume model by default, but it can be
// changed into a DG formulation by changing the fe_order and
// fe_family parameters in the input file appropriately. The example
// uses the Lax-Friedrichs numerical flux. This is known to be more
// diffusive than e.g. upwinding, but it is also relatively simple
// to implement since one does not need to consider the advective
// velocity direction during assembly. Finally, the example comes
// with two different time discretization options: explicit "forward"
// Euler and explicit fourth-order Runge-Kutta (RK4). The latter
// is used by default, but the "temporal_discretization_type" can
// be changed to ForwardEuler in the input file to test that option.

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"

// Application includes
#include "advection_system.h"

// C++ include files that we need
#include <iostream>

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function to set the initial condition
void set_initial_condition (EquationSystems& es,
                            const std::string& system_name);
Number initial_condition (const Point& p,
                          const Parameters& parameters,
                          const std::string&,
                          const std::string&);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Create GetPot object to parse the command line
  GetPot command_line (argc, argv);

  // Allow the user specify a custom input file on the command line
  // using the "-i" flag, and provide a sensible default value in case
  // they don't.
  std::string parameters_filename = "advection_2D.in";
  if ( command_line.search(2, "-i", "--input_file") )
    parameters_filename = command_line.next(parameters_filename);

  // Read mesh options from the input file
  GetPot infile(parameters_filename);
  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions

  // Build a mesh
  Mesh mesh (init.comm());
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       0., 1.,
                                       0., 1.,
                                       QUAD9);

  // Create an equation systems object
  EquationSystems equation_systems (mesh);

  // Create an AdvectionSystem within the EquationSystems
  AdvectionSystem & advection_system =
    equation_systems.add_system<AdvectionSystem> ("Advection2D");

  // Process the input file values
  advection_system.process_parameters_file(parameters_filename);

  // Set the initial condition function
  advection_system.attach_init_function (set_initial_condition);

  equation_systems.init ();

  equation_systems.print_info();
  mesh.print_info();

  advection_system.print_info();

  advection_system.assemble_all_matrices();
  // advection_system.write_out_discretization_matrices();

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems ("claw_solution.0000.e", equation_systems);
#endif

  advection_system.solve_conservation_law();

  return 0;
}

void set_initial_condition (EquationSystems& es,
                            const std::string&)
{
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem>("Advection2D");

  system.project_solution(initial_condition, /*gradient=*/nullptr, es.parameters);
}

Number initial_condition (const Point& p,
                          const Parameters& ,
                          const std::string&,
                          const std::string&)
{
  Real x = p(0);
  Real y = p(1);

  return std::exp( -100.*(std::pow(x-0.5,2.) + std::pow(y-0.5,2.)) );
}
