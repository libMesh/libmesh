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


// <h1>FEMSystem Example 2</h1>
// \author Robert Weidlich
// \date 2012

#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/time_solver.h"
#include "libmesh/transient_system.h"
#include "libmesh/vtk_io.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

using namespace libMesh;

#include "solid_system.h"

void setup(EquationSystems & systems,
           Mesh & mesh,
           GetPot & args)
{
  const unsigned int dim = mesh.mesh_dimension();
  // We currently invert tensors with the assumption that they're 3x3
  libmesh_assert (dim == 3);

  // Generating Mesh
  ElemType eltype = Utility::string_to_enum<ElemType>(args("mesh/generation/element_type", "hex8"));
  int nx = args("mesh/generation/num_elem", 4, 0);
  int ny = args("mesh/generation/num_elem", 4, 1);
  int nz = dim > 2 ? args("mesh/generation/num_elem", 4, 2) : 0;
  double origx = args("mesh/generation/origin", -1.0, 0);
  double origy = args("mesh/generation/origin", -1.0, 1);
  double origz = args("mesh/generation/origin", 0.0, 2);
  double sizex = args("mesh/generation/size", 2.0, 0);
  double sizey = args("mesh/generation/size", 2.0, 1);
  double sizez = args("mesh/generation/size", 2.0, 2);
  MeshTools::Generation::build_cube(mesh, nx, ny, nz,
                                    origx, origx+sizex, origy, origy+sizey, origz, origz+sizez, eltype);

  // Creating Systems
  SolidSystem & imms = systems.add_system<SolidSystem> ("solid");
  imms.args = args;

  // Build up auxiliary system
  ExplicitSystem & aux_sys = systems.add_system<TransientExplicitSystem>("auxiliary");

  // Initialize the system
  systems.parameters.set<unsigned int>("phase") = 0;
  systems.init();

  imms.save_initial_mesh();

  // Fill global solution vector from local ones
  aux_sys.current_local_solution->close();
  *aux_sys.solution = *aux_sys.current_local_solution;
  aux_sys.solution->close();
}



void run_timestepping(EquationSystems & systems, GetPot & args)
{
  TransientExplicitSystem & aux_system = systems.get_system<TransientExplicitSystem>("auxiliary");

  SolidSystem & solid_system = systems.get_system<SolidSystem>("solid");

  std::unique_ptr<VTKIO> io = libmesh_make_unique<VTKIO>(systems.get_mesh());

  Real duration = args("duration", 1.0);

  for (unsigned int t_step = 0; t_step < duration/solid_system.deltat; t_step++)
    {
      // Progress in current phase [0..1]
      Real progress = t_step * solid_system.deltat / duration;
      systems.parameters.set<Real>("progress") = progress;
      systems.parameters.set<unsigned int>("step") = t_step;

      // Update message

      out << "===== Time Step " << std::setw(4) << t_step;
      out << " (" << std::fixed << std::setprecision(2) << std::setw(6) << (progress*100.) << "%)";
      out << ", time = " << std::setw(7) << solid_system.time;
      out << " =====" << std::endl;

      // Advance timestep in auxiliary system
      aux_system.current_local_solution->close();
      aux_system.old_local_solution->close();
      *aux_system.older_local_solution = *aux_system.old_local_solution;
      *aux_system.old_local_solution = *aux_system.current_local_solution;

      out << "Solving Solid" << std::endl;
      solid_system.solve();
      aux_system.current_local_solution->close();
      *aux_system.solution = *aux_system.current_local_solution;
      aux_system.solution->close();

      // Carry out the adaptive mesh refinement/coarsening
      out << "Doing a reinit of the equation systems" << std::endl;
      systems.reinit();

      if (t_step % args("output/frequency", 1) == 0)
        {
          std::stringstream file_name;
          file_name << args("results_directory", "./")
                    << "fem_"
                    << std::setw(6)
                    << std::setfill('0')
                    << t_step
                    << ".pvtu";

          io->write_equation_systems(file_name.str(), systems);
        }
      // Advance to the next timestep in a transient problem
      out << "Advancing to next step" << std::endl;
      solid_system.time_solver->advance_timestep();
    }
}



int main(int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip this example if we do not meet certain requirements
#ifndef LIBMESH_HAVE_VTK
  libmesh_example_requires(false, "--enable-vtk");
#endif

  // Trilinos gives us an inverted element on this one...
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc");

  // Threaded assembly doesn't currently work with the moving mesh
  // code.
  // We'll skip this example for now.
  if (libMesh::n_threads() > 1)
    {
      libMesh::out << "We skip fem_system_ex2 when using threads.\n"
                   << std::endl;
      return 0;
    }

  // read simulation parameters from file
  GetPot args = GetPot("solid.in");

  // Create System and Mesh
  int dim = args("mesh/generation/dimension", 3);
  libmesh_example_requires(dim <= LIBMESH_DIM, "3D support");

  // Create a mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);

  EquationSystems systems(mesh);

  // Create and set systems up
  setup(systems, mesh, args);

  // run the systems
  run_timestepping(systems, args);

  out << "Finished calculations" << std::endl;
  return 0;
}
