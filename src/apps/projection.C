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


// Open the input mesh and corresponding solution file named in command line
// arguments, open the output mesh, project that solution onto the
// output mesh, and write a corresponding output solution file.

#include <map>
#include <string>

#include "libmesh/libmesh.h"

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"
#include "libmesh/replicated_mesh.h"


using namespace libMesh;


// If there's a missing input argument, then print a help message
void usage_error(const char * progname)
{
  libMesh::out << "Options: " << progname << '\n'
               << " --dim d               mesh dimension           [default: autodetect]\n"
               << " --inmesh filename     input mesh file\n"
               << " --insoln filename     input solution file\n"
               << " --outmesh filename    output mesh file         [default: out_<inmesh>]\n"
               << " --outsoln filename    output solution file     [default: out_<insoln>]\n"
               << std::endl;

  exit(1);
}

// Get an input argument, or print a help message if it's missing
template <typename T>
T assert_argument (GetPot & cl,
                   const std::string & argname,
                   const char * progname,
                   const T & defaultarg)
{
  if (!cl.search(argname))
    {
      libMesh::err << ("No " + argname + " argument found!") << std::endl;
      usage_error(progname);
    }
  return cl.next(defaultarg);
}

// Global collections of MeshFunctions are necessary to work with
// global functions fptr and gptr.
// TODO: libMesh needs functor-based alternatives to these types of
// function arguments
std::string current_sys_name;
std::map<std::string, MeshFunction *> mesh_functions;

// Return the function value on the old mesh and solution
Number fptr(const Point & p,
            const Parameters &,
            const std::string & libmesh_dbg_var(sys_name),
            const std::string & unknown_name)
{
  libmesh_assert_equal_to (sys_name, current_sys_name);
  libmesh_assert(mesh_functions.count(unknown_name));
  libmesh_assert(mesh_functions[unknown_name]);

  MeshFunction & meshfunc = *mesh_functions[unknown_name];

  return meshfunc(p);
}

// Return the function gradient on the old mesh and solution
Gradient gptr(const Point & p,
              const Parameters &,
              const std::string & libmesh_dbg_var(sys_name),
              const std::string & unknown_name)
{
  libmesh_assert_equal_to (sys_name, current_sys_name);
  libmesh_assert(mesh_functions.count(unknown_name));
  libmesh_assert(mesh_functions[unknown_name]);

  MeshFunction & meshfunc = *mesh_functions[unknown_name];

  return meshfunc.gradient(p);
}


int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  GetPot cl(argc, argv);

  // In case the mesh file doesn't let us auto-infer dimension, we let
  // the user specify it on the command line
  const unsigned char requested_dim =
    cast_int<unsigned char>(cl.follow(3, "--dim"));

  // Load the old mesh from --inmesh filename.
  // Keep it serialized; we don't want elements on the new mesh to be
  // looking for data on old mesh elements that live off-processor.
  ReplicatedMesh old_mesh(init.comm(), requested_dim);

  const std::string meshname =
    assert_argument(cl, "--inmesh", argv[0], std::string("mesh.xda"));

  old_mesh.read(meshname);
  std::cout << "Old Mesh:" << std::endl;
  old_mesh.print_info();

  // Load the new mesh from --outmesh filename
  Mesh new_mesh(init.comm(), requested_dim);

  const std::string outmeshname = cl.follow(std::string("out_"+meshname), "--outmesh");

  new_mesh.read(outmeshname);
  std::cout << "New Mesh:" << std::endl;
  new_mesh.print_info();

  // Load the old solution from --insoln filename
  // Construct the new solution from the old solution's headers, so
  // that we get the system names and types, variable names and types,
  // etc.
  const std::string solnname =
    assert_argument(cl, "--insoln", argv[0], std::string("soln.xda"));

  EquationSystems old_es(old_mesh);
  EquationSystems new_es(new_mesh);

  XdrMODE read_mode;

  if (solnname.rfind(".xdr") < solnname.size())
    read_mode = DECODE;
  else if (solnname.rfind(".xda") < solnname.size())
    read_mode = READ;
  else
    libmesh_error_msg("Unrecognized file extension on " << solnname);

  old_es.read(solnname, read_mode,
              EquationSystems::READ_HEADER |
              EquationSystems::READ_DATA |
              EquationSystems::READ_ADDITIONAL_DATA |
              EquationSystems::READ_BASIC_ONLY);

  new_es.read(solnname, read_mode,
              EquationSystems::READ_HEADER |
              EquationSystems::READ_BASIC_ONLY);

  old_es.print_info();

  unsigned int n_systems = old_es.n_systems();
  libmesh_assert_equal_to (new_es.n_systems(), n_systems);

  // For each system, serialize the solution so we can project it onto
  // a potentially-very-different partitioning
  for (unsigned int i = 0; i != n_systems; ++i)
    {
      System & old_sys = old_es.get_system(i);
      current_sys_name = old_sys.name();

      libMesh::out << "Projecting system " << current_sys_name << std::endl;

      libmesh_assert (new_es.has_system(current_sys_name));

      System & new_sys = new_es.get_system(current_sys_name);
      unsigned int n_vars = old_sys.n_vars();
      libmesh_assert_equal_to (new_sys.n_vars(), n_vars);

      std::unique_ptr<NumericVector<Number>> comparison_soln =
        NumericVector<Number>::build(old_sys.comm());
      std::vector<Number> global_soln;
      old_sys.update_global_solution(global_soln);
      comparison_soln->init(old_sys.solution->size(), true, SERIAL);
      (*comparison_soln) = global_soln;

      // For each variable, construct a MeshFunction returning that
      // variable's value
      for (unsigned int j = 0; j != n_vars; ++j)
        {
          libMesh::out << " with variable " << old_sys.variable_name(j) << std::endl;

          MeshFunction * mesh_func =
            new MeshFunction(old_es, *comparison_soln,
                             old_sys.get_dof_map(), j);
          mesh_func->init();
          mesh_functions[old_sys.variable_name(j)] = mesh_func;
        }

      // Project all variables to the new system
      new_sys.project_solution(fptr, gptr, old_es.parameters);

      // Clean up the MeshFunctions here so we don't bloat memory
      for (unsigned int j = 0; j != n_vars; ++j)
        delete mesh_functions[old_sys.variable_name(j)];
    }

  // Write out the new solution file
  const std::string outsolnname = cl.follow(std::string("out_"+solnname), "--outsoln");

  new_es.write(outsolnname.c_str(),
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
  libMesh::out << "Wrote solution " << outsolnname << std::endl;

  return 0;
}
