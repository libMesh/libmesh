// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Solution Transfer Example 1 - </h1>
// \author Derek Gaston
// \date 2013
//
// This example demonstrates how to use the DTKSolutionTransfer object
// to transfer a solution from one Mesh to another.

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK
#include "libmesh/dtk_solution_transfer.h"
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;


Number initial_value(const Point & p,
                     const Parameters & /* parameters */,
                     const std::string &,
                     const std::string &)
{
  return p(0)*p(0) + 1; // x^2 + 1
}

void initialize(EquationSystems & es,
                const std::string & system_name)
{
  ExplicitSystem & system = es.get_system<ExplicitSystem>(system_name);
  es.parameters.set<Real> ("time") = system.time = 0;
  system.project_solution(initial_value, libmesh_nullptr, es.parameters);
}

int main(int argc, char * argv[])
{
  LibMeshInit init (argc, argv);

#if !defined(LIBMESH_TRILINOS_HAVE_DTK)
  // Skip this example (and use a different return code) if libMesh
  // was compiled without Trilinos+DTK support.
  libmesh_example_requires(false, "--enable-trilinos");

#else

  Mesh from_mesh(init.comm());
  MeshTools::Generation::build_cube(from_mesh, 4, 4, 4, 0, 1, 0, 1, 0, 1, HEX8);
  from_mesh.print_info();
  EquationSystems from_es(from_mesh);
  System & from_sys = from_es.add_system<ExplicitSystem>("From");
  unsigned int from_var = from_sys.add_variable("from");
  from_sys.attach_init_function(initialize);
  from_es.init();

  ExodusII_IO(from_mesh).write_equation_systems("from.e", from_es);

  Mesh to_mesh(init.comm());
  MeshTools::Generation::build_cube(to_mesh, 5, 5, 5, 0, 1, 0, 1, 0, 1, TET4);
  to_mesh.print_info();
  EquationSystems to_es(to_mesh);
  System & to_sys = to_es.add_system<ExplicitSystem>("To");
  unsigned int to_var = to_sys.add_variable("to");
  to_es.init();

  DTKSolutionTransfer dtk_transfer(init.comm());

  dtk_transfer.transfer(from_sys.variable(from_var), to_sys.variable(to_var));

  to_es.update();
  ExodusII_IO(to_mesh).write_equation_systems("to.e", to_es);

#endif

  return 0;
}
