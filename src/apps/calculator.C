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


// Open the input mesh and corresponding solution file named in command line
// arguments, parse the function specified in a command line argument,
// L2-project its value onto the mesh, and output the new solution.


// libMesh includes
#include "L2system.h"
#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/point.h"
#include "libmesh/steady_solver.h"
#include "libmesh/wrapped_functor.h"
#include "libmesh/elem.h"
#include "libmesh/parallel_implementation.h"
#include "libmesh/string_to_enum.h"


// C++ includes
#include <memory>
#include <string>


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
               << " --subdomain num       subdomain id restriction [default: all subdomains]\n"
               << " --integral            only calculate integral, not projection\n"
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


struct Integrate {
  Integrate (const System & sys,
             const FEMFunctionBase<Number> & f) :
    _sys(sys), _f(f.clone()), _integral(0) {}

  Integrate (Integrate & other, Threads::split) :
    _sys(other._sys), _f(other._f->clone()), _integral(0) {}

  void operator() (const ConstElemRange & range) {
    FEMContext context(_sys);

    FEBase * elem_fe = nullptr;
    context.get_element_fe(0, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();
    const std::vector<Point> & xyz = elem_fe->get_xyz();

    _f->init_context(context);

    // If _f didn't need anything, that wasn't a fluke
    elem_fe->get_nothing();

    // We only integrate on the highest dimensional elements in the
    // mesh, lest we be adding apples to oranges
    const unsigned int mesh_dim = _sys.get_mesh().mesh_dimension();

    for (const auto & elem : range)
      {
        if (elem->dim() < mesh_dim)
          continue;

        context.pre_fe_reinit(_sys, elem);
        context.elem_fe_reinit();

        for (auto qp : index_range(JxW))
          {
            const Number output = (*_f)(context, xyz[qp]);

            _integral += output * JxW[qp];
          }
      }
  }

  void join (const Integrate & other)
  {
    _integral += other._integral;
  }

  Number integral () const { return _integral; }

private:
  const System & _sys;

  std::unique_ptr<FEMFunctionBase<Number>> _f;

  Number _integral;
};


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
    assert_argument(cl, "--inmesh", argv[0], std::string(""));

  const std::string solnname = cl.follow(std::string(""), "--insoln");

  if (meshname != "")
    {
      old_mesh.read(meshname);
      std::cout << "Old Mesh:" << std::endl;
      old_mesh.print_info();
    }

  // Load the old solution from --insoln filename, if that's been
  // specified.
  EquationSystems old_es(old_mesh);
  std::string current_sys_name = "new_sys";

  const std::string calcfunc =
    assert_argument(cl, "--calc", argv[0], std::string(""));

  const std::string family =
    cl.follow(std::string("LAGRANGE"), "--family");

  const int order = cl.follow(1, "--order");

  std::unique_ptr<FEMFunctionBase<Number>> goal_function;

  if (solnname != "")
    {
      old_es.read(solnname,
                  EquationSystems::READ_HEADER |
                  EquationSystems::READ_DATA |
                  EquationSystems::READ_ADDITIONAL_DATA |
                  EquationSystems::READ_BASIC_ONLY);

      old_es.print_info();

      const unsigned int sysnum =
        cl.follow(0, "--insys");

      libmesh_assert_less(sysnum, old_es.n_systems());

      System & old_sys = old_es.get_system(sysnum);

      current_sys_name = old_sys.name();

      goal_function =
        std::make_unique<ParsedFEMFunction<Number>>(old_sys, calcfunc);
    }
  else
    {
      current_sys_name = "bare_sys";

      // FEMContext doesn't like seeing systems with no variables
      System & dummy_sys = old_es.add_system<System>(current_sys_name);

      // And if we're using an IsoGeometric Analysis mesh, we'd better
      // make sure the user can specify order and family to make it
      // iso.
      dummy_sys.add_variable("dummy", static_cast<Order>(order),
                             Utility::string_to_enum<FEFamily>(family));

      old_es.init();

      goal_function =
        std::make_unique<WrappedFunctor<Number>>(ParsedFunction<Number>(calcfunc));
    }

  libMesh::out << "Calculating with system " << current_sys_name << std::endl;

  // Subdomains to integrate on
  std::set<libMesh::subdomain_id_type> subdomains_list;
  cl.disable_loop();
  while (cl.search(1, "--subdomain"))
    {
      subdomain_id_type tmp = Elem::invalid_subdomain_id;
      tmp = cl.next(tmp);
      subdomains_list.insert(tmp);
    }
  cl.enable_loop();

  std::string default_outsolnname = "out_soln.xda";
  if (solnname != "")
    default_outsolnname = "out_"+solnname;
  const std::string outsolnname =
    cl.follow(default_outsolnname, "--outsoln");

  if (!cl.search("--integral"))
    {
      // Create a new mesh and a new EquationSystems
      Mesh new_mesh(init.comm(), requested_dim);
      new_mesh.read(meshname);

      EquationSystems new_es(new_mesh);

      L2System & new_sys = new_es.add_system<L2System>(current_sys_name);

      new_sys.subdomains_list() = std::move(subdomains_list);

      new_sys.time_solver =
        std::make_unique<SteadySolver>(new_sys);

      new_sys.fe_family() = family;
      new_sys.fe_order() = order;

      new_sys.goal_func = goal_function->clone();

      if (solnname != "")
        new_sys.input_system = &old_es.get_system(current_sys_name);

      new_es.init();

      DiffSolver & solver = *(new_sys.time_solver->diff_solver().get());
      solver.quiet = false;
      solver.verbose = true;
      solver.relative_step_tolerance = 1e-10;

      new_sys.solve();

      // Write out the new solution file
      new_es.write(outsolnname.c_str(),
                   EquationSystems::WRITE_DATA |
                   EquationSystems::WRITE_ADDITIONAL_DATA);
      libMesh::out << "Wrote solution " << outsolnname << std::endl;
    }
  else
    {
      // If we aren't doing a projection, then we just do an integral
      System & old_sys = old_es.get_system<System>(current_sys_name);

      ConstElemRange active_local_elem_range
        (old_mesh.active_local_elements_begin(),
         old_mesh.active_local_elements_end());

      Integrate integrate(old_sys, *goal_function);

      Threads::parallel_reduce (active_local_elem_range,
                                integrate);

      Number integral = integrate.integral();
      old_mesh.comm().sum(integral);
      std::cout << "Integral is " << integral << std::endl;
    }

  return 0;
}
