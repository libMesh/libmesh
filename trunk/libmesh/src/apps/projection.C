// Open the input mesh and corresponding solution file named in command line
// arguments, open the output mesh, project that solution onto the
// output mesh, and write a corresponding output solution file.

#include <map>
#include <string>

#include "libmesh.h"

#include "dof_map.h"
#include "equation_systems.h"
#include "getpot.h"
#include "mesh.h"
#include "mesh_function.h"
#include "numeric_vector.h"
#include "point.h"


// If there's a missing input argument, then print a help message
void usage_error(const char *progname)
{
  std::cout << "Options: " << progname << '\n'
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
T assert_argument (GetPot &cl,
                   const std::string &argname,
                   const char        *progname,
                   const T&          defaultarg)
{
  if(!cl.search(argname))
    {
      std::cerr << ("No " + argname + " argument found!") << std::endl;
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
Number fptr(const Point& p,
            const Parameters&,
            const std::string& sys_name,
            const std::string& unknown_name)
{
  libmesh_assert(sys_name == current_sys_name);
  libmesh_assert(mesh_functions.count(unknown_name));
  libmesh_assert(mesh_functions[unknown_name]);

  MeshFunction &meshfunc = *mesh_functions[unknown_name];

  return meshfunc(p);
}

// Return the function gradient on the old mesh and solution
Gradient gptr(const Point& p,
              const Parameters&,
              const std::string& sys_name,
              const std::string& unknown_name)
{
  libmesh_assert(sys_name == current_sys_name);
  libmesh_assert(mesh_functions.count(unknown_name));
  libmesh_assert(mesh_functions[unknown_name]);

  MeshFunction &meshfunc = *mesh_functions[unknown_name];

  return meshfunc.gradient(p);
}


int main(int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  GetPot cl(argc, argv);

  // In case the mesh file doesn't let us auto-infer dimension, we let
  // the user specify it on the command line
  const unsigned int requested_dim = cl.follow(3, "--dim");

  // Load the old mesh from --inmesh filename
  Mesh old_mesh(requested_dim);

  const std::string meshname =
    assert_argument(cl, "--inmesh", argv[0], std::string("mesh.xda"));

  old_mesh.read(meshname);
  old_mesh.print_info();

  // Load the new mesh from --outmesh filename
  Mesh new_mesh(requested_dim);

  const std::string outmeshname = cl.follow(std::string("out_"+meshname), "--outmesh");

  new_mesh.read(outmeshname);

  // Load the old solution from --insoln filename
  // Construct the new solution from the old solution's headers, so
  // that we get the system names and types, variable names and types,
  // etc.
  const std::string solnname =
    assert_argument(cl, "--insoln", argv[0], std::string("soln.xda"));

  EquationSystems old_es(old_mesh);
  EquationSystems new_es(new_mesh);

  libMeshEnums::XdrMODE read_mode;

  if (solnname.rfind(".xdr") < solnname.size())
    read_mode = DECODE;
  else if (solnname.rfind(".xda") < solnname.size())
    read_mode = READ;
  else
    {
      std::cerr << "Unrecognized file extension on " << solnname << std::endl;
      libmesh_error();
    }

  old_es.read(solnname, read_mode,
              EquationSystems::READ_HEADER |
              EquationSystems::READ_DATA |
              EquationSystems::READ_ADDITIONAL_DATA);

  new_es.read(solnname, read_mode,
              EquationSystems::READ_HEADER);

  old_es.print_info();

  unsigned int n_systems = old_es.n_systems();
  libmesh_assert(new_es.n_systems() == n_systems);

  // For each system, serialize the solution so we can project it onto
  // a potentially-very-different partitioning
  for (unsigned int i = 0; i != n_systems; ++i)
    {
      System &old_sys = old_es.get_system(i);
      current_sys_name = old_sys.name();

      libmesh_assert (new_es.has_system(current_sys_name));

      System &new_sys = new_es.get_system(current_sys_name);
      unsigned int n_vars = old_sys.n_vars();
      libmesh_assert(new_sys.n_vars() == n_vars);

      AutoPtr<NumericVector<Number> > comparison_soln =
        NumericVector<Number>::build();
      std::vector<Number> global_soln;
      old_sys.update_global_solution(global_soln);
      comparison_soln->init(old_sys.solution->size(), true, SERIAL);
      (*comparison_soln) = global_soln;

      // For each variable, construct a MeshFunction returning that
      // variable's value
      for (unsigned int j = 0; j != n_vars; ++j)
        {
          MeshFunction *mesh_func =
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
  std::cout << "Wrote solution " << outsolnname << std::endl;

  return 0;
}
