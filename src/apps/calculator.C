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


// Open the input mesh and corresponding solution file named in command line
// arguments, parse the function specified in a command line argument,
// L2-project its value onto the mesh, and output the new solution.

#include <map>
#include <string>

#include "L2system.h"

#include "libmesh/libmesh.h"

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/point.h"
#include "libmesh/steady_solver.h"


using namespace libMesh;


// If there's a missing input argument, then print a help message
void usage_error(const char * progname)
{
  libMesh::out << "Options: " << progname << '\n'
               << " --dim d               mesh dimension           [default: autodetect]\n"
               << " --inmesh    filename  input mesh file\n"
               << " --insoln    filename  input solution file\n"
               << " --calc      func      function to calculate\n"
               << " --insys     num       input system number      [default: 0]\n"
               << " --outsoln   filename  output solution file     [default: out_<insoln>]\n"
               << " --family    famname   output FEM family        [default: LAGRANGE]\n"
               << " --order     num       output FEM order         [default: 1]\n"
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
  Mesh old_mesh(init.comm(), requested_dim);

  const std::string meshname =
    assert_argument(cl, "--inmesh", argv[0], std::string("mesh.xda"));

  old_mesh.read(meshname);
  std::cout << "Mesh:" << std::endl;
  old_mesh.print_info();

  // Create a new mesh for a new EquationSystems
  Mesh new_mesh(init.comm(), requested_dim);
  new_mesh.read(meshname);

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

  old_es.print_info();


  const unsigned int sysnum =
    cl.follow(0, "--insys");

  libmesh_assert_less(sysnum, old_es.n_systems());

  System & old_sys = old_es.get_system(sysnum);
  std::string current_sys_name = old_sys.name();

  libMesh::out << "Calculating with system " << current_sys_name << std::endl;

  L2System & new_sys = new_es.add_system<L2System>(current_sys_name);

  new_sys.time_solver =
    libmesh_make_unique<SteadySolver>(new_sys);

  new_sys.fe_family() =
    cl.follow(std::string("LAGRANGE"), "--family");

  new_sys.fe_order() =
    cl.follow(1, "--order");

  const std::string calcfunc =
    assert_argument(cl, "--calc", argv[0], std::string(""));

  ParsedFEMFunction<Number> goal_function(old_sys, calcfunc);

  new_sys.goal_func = goal_function.clone();
  new_sys.input_system = &old_sys;

  new_es.init();

  DiffSolver & solver = *(new_sys.time_solver->diff_solver().get());
  solver.quiet = false;
  solver.verbose = true;
  solver.relative_step_tolerance = 1e-10;

  new_sys.solve();

  // Write out the new solution file
  const std::string outsolnname = cl.follow(std::string("out_"+solnname), "--outsoln");

  new_es.write(outsolnname.c_str(),
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
  libMesh::out << "Wrote solution " << outsolnname << std::endl;

  return 0;
}
