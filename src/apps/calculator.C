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


// Open the input mesh and corresponding solution file named in command line
// arguments, parse the function specified in a command line argument,
// Project its value onto the mesh, and output the new solution.


// libMesh includes
#include "L2system.h"
#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/point.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/steady_solver.h"
#include "libmesh/wrapped_functor.h"
#include "libmesh/elem.h"
#include "libmesh/parallel_implementation.h"
#include "libmesh/string_to_enum.h"

// For error integration
#include "libmesh/error_vector.h"
#include "libmesh/exact_solution.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/fourth_error_estimators.h"

// C++ includes
#include <memory>
#include <string>


using namespace libMesh;


// If there's a missing input argument, then print a help message
void usage_error(const char * progname)
{
  libMesh::out << "Options: " << progname << '\n'
               << " --dim        d         mesh dimension               [default: autodetect]\n"
               << " --inmesh     filename  input mesh file\n"
               << " --inmat      filename  input constraint matrix      [default: none]\n"
               << " --mattol     filename  constraint tolerance when testing mesh connectivity\n"
               << "                                                     [default: 0]\n"
               << " --insoln     filename  input solution file\n"
               << " --calc       func      function to calculate\n"
               << " --insys      sysnum    input system number          [default: 0]\n"
               << " --outsoln    filename  output solution file         [default: out_<insoln>]\n"
               << " --family     famname   output FEM family            [default: LAGRANGE]\n"
               << " --order      p         output FEM order             [default: 1]\n"
               << " --subdomain  \"sbd_ids\" each subdomain to check      [default: all subdomains]\n"
               << " --hilbert    order     Hilbert space to project in  [default: 0 => H0]\n"
               << " --fdm_eps    eps       Central diff for dfunc/dx    [default: " << TOLERANCE << "]\n"
               << " --error_q    extra_q   integrate projection error, with adjusted\n"
               << "                       (extra) quadrature order      [default: off, suggested: 0]\n"
               << " --jump_error order     calculate jump error indicator, for specified\n"
               << "                       Hilbert order                 [default: off]\n"
               << " --jump_slits           calculate jumps across slits [default: off]\n"
               << " --integral             only calculate func integral, not projection\n"
               << std::endl;

  exit(1);
}

// Get an input argument, or print a help message if it's missing
template <typename T>
T assert_argument (const std::string & argname,
                   const char * progname,
                   const T & defaultarg)
{
  if (!libMesh::on_command_line(argname))
    {
      libMesh::err << ("No " + argname + " argument found!") << std::endl;
      usage_error(progname);
    }
  return libMesh::command_line_next(argname, defaultarg);
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

  // In case the mesh file doesn't let us auto-infer dimension, we let
  // the user specify it on the command line
  const unsigned char requested_dim =
    cast_int<unsigned char>(libMesh::command_line_next("--dim", 3));

  // Load the old mesh from --inmesh filename.
  // Keep it serialized; we don't want elements on the new mesh to be
  // looking for data on old mesh elements that live off-processor.
  Mesh old_mesh(init.comm(), requested_dim);

  const std::string meshname =
    assert_argument("--inmesh", argv[0], std::string(""));

  libMesh::out << "Reading mesh " << meshname << std::endl;
  old_mesh.read(meshname);

  const std::string matname =
    libMesh::command_line_next("--inmat", std::string(""));

  if (matname != "")
    {
      libMesh::out << "Reading matrix " << matname << std::endl;

      // For extraction matrices Coreform has been experimenting with
      // PETSc solvers which take the transpose of what we expect, so
      // we'll un-transpose here.
      auto matrix = SparseMatrix<Number>::build (old_mesh.comm());
      matrix->read(matname);
      matrix->get_transpose(*matrix);

      old_mesh.copy_constraint_rows(*matrix);
    }

  libMesh::out << "Mesh:" << std::endl;
  old_mesh.print_info();

  // If we're not using a distributed mesh, this is cheap info to add
  if (old_mesh.is_serial_on_zero())
    {
      const Real mat_tol =
        libMesh::command_line_next("--mattol", Real(0));

      const dof_id_type n_components =
        MeshTools::n_connected_components(old_mesh, mat_tol);
      libMesh::out << "Mesh has " << n_components << " connected components." << std::endl;
    }

  const std::string solnname =
    libMesh::command_line_next("--insoln", std::string(""));

  // Load the old solution from --insoln filename, if that's been
  // specified.
  EquationSystems old_es(old_mesh);
  std::string current_sys_name = "new_sys";

  const std::string calcfunc =
    assert_argument("--calc", argv[0], std::string(""));

  const std::string family =
    libMesh::command_line_next("--family", std::string("LAGRANGE"));

  const unsigned int order = libMesh::command_line_next("--order", 1u);

  std::unique_ptr<FEMFunctionBase<Number>> goal_function;

  if (solnname != "")
    {
      libMesh::out << "Reading solution " << solnname << std::endl;

      old_es.read(solnname,
                  EquationSystems::READ_HEADER |
                  EquationSystems::READ_DATA |
                  EquationSystems::READ_ADDITIONAL_DATA |
                  EquationSystems::READ_BASIC_ONLY);

      old_es.print_info();

      const unsigned int sysnum =
        libMesh::command_line_next("--insys", 0);

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

      old_es.print_info();

      goal_function =
        std::make_unique<WrappedFunctor<Number>>(ParsedFunction<Number>(calcfunc));
    }

  libMesh::out << "Calculating with system " << current_sys_name << std::endl;

  // Subdomains to integrate on
  std::vector<libMesh::subdomain_id_type> subdomain_vec;
  libMesh::command_line_vector("--subdomain", subdomain_vec);
  std::set<libMesh::subdomain_id_type> subdomains_list(subdomain_vec.begin(),
                                                       subdomain_vec.end());

  std::string default_outsolnname = "out_soln.xda";
  if (solnname != "")
    default_outsolnname = "out_"+solnname;
  const std::string outsolnname =
    libMesh::command_line_next("--outsoln", default_outsolnname);

  // Output results in high precision
  libMesh::out << std::setprecision(std::numeric_limits<Real>::max_digits10);

  if (!libMesh::on_command_line("--integral"))
    {
      // Create a new mesh and a new EquationSystems
      Mesh new_mesh(init.comm(), requested_dim);
      new_mesh.read(meshname);

      EquationSystems new_es(new_mesh);

      HilbertSystem & new_sys = new_es.add_system<HilbertSystem>(current_sys_name);

      new_sys.hilbert_order() = libMesh::command_line_next("--hilbert", 0u);

      new_sys.subdomains_list() = std::move(subdomains_list);

      new_sys.time_solver =
        std::make_unique<SteadySolver>(new_sys);

      new_sys.fe_family() = family;
      new_sys.fe_order() = order;

      new_sys.set_goal_func(*goal_function);

      const Real fdm_eps = libMesh::command_line_next("--fdm_eps", Real(TOLERANCE));

      new_sys.set_fdm_eps(fdm_eps);

      if (solnname != "")
        new_sys.input_system = &old_es.get_system(current_sys_name);

      new_es.init();

      DiffSolver & solver = *(new_sys.time_solver->diff_solver().get());
      solver.quiet = false;
      solver.verbose = true;
      solver.relative_step_tolerance = 1e-10;

      new_sys.solve();

      // Integrate the error if requested
      if (libMesh::on_command_line("--error_q"))
        {
          const unsigned int error_q =
            libMesh::command_line_next("--error_q", 0u);

          // We just add "u" now but maybe we'll change that
          for (auto v : make_range(new_sys.n_vars()))
            {
              ExactSolution exact_sol(new_es);
              exact_sol.attach_exact_value(0, goal_function.get());
              FDMGradient<Gradient> fdm_gradient(*goal_function, fdm_eps);
              exact_sol.attach_exact_deriv(0, &fdm_gradient);
              exact_sol.extra_quadrature_order(error_q);

              const std::string var_name = new_sys.variable_name(v);
              exact_sol.compute_error(current_sys_name, var_name);

              libMesh::out << "L2 norm of " << var_name << ": " <<
                new_sys.calculate_norm(*new_sys.solution, v, L2) <<
                std::endl;

              libMesh::out << "L2 error in " << var_name << ": " <<
                  exact_sol.l2_error(current_sys_name, var_name) <<
                  std::endl;

              if (new_sys.hilbert_order() > 0)
                {
                  libMesh::out << "H1 norm of " << var_name << ": " <<
                    new_sys.calculate_norm(*new_sys.solution, v, H1) <<
                    std::endl;

                  libMesh::out << "L2 error in " << var_name << ": " <<
                      exact_sol.h1_error(current_sys_name, var_name) <<
                      std::endl;
                }
            }
        }

      // Calculate the jump error indicator if requested
      if (libMesh::on_command_line("--jump_error"))
        {
          const unsigned int jump_error_hilbert =
            libMesh::command_line_next("--jump_error", 0u);

          for (auto v : make_range(new_sys.n_vars()))
            {
              std::unique_ptr<JumpErrorEstimator> error_estimator;
              SystemNorm error_norm;
              error_norm.set_weight(0,0);
              error_norm.set_weight(v+1,0);
              error_norm.set_weight(v,1.0);

              if (jump_error_hilbert == 0)
                {
                  error_estimator = std::make_unique<DiscontinuityMeasure>();
                  error_norm.set_type(v, L2);
                }
              else if (jump_error_hilbert == 1)
                {
                  error_estimator = std::make_unique<KellyErrorEstimator>();
                  error_norm.set_type(v, H1_SEMINORM);
                }
              else if (jump_error_hilbert == 2)
                {
                  error_estimator = std::make_unique<LaplacianErrorEstimator>();
                  error_norm.set_type(v, H2_SEMINORM);
                }
              else
                libmesh_not_implemented();

              if (libMesh::on_command_line("--jump_slits"))
                error_estimator->integrate_slits = true;

              ErrorVector error_per_cell;
              error_estimator->error_norm = error_norm;
              error_estimator->estimate_error(new_sys, error_per_cell);

              const std::string var_name = new_sys.variable_name(v);

              // I really want to just stop supporting libstdc++-8.
              // const Real total_error = std::reduce(error_per_cell.begin(), error_per_cell.end());

              ErrorVectorReal total_error = 0;
              for (auto cell_error : error_per_cell)
                total_error += cell_error;

              libMesh::out << "H" << jump_error_hilbert << " error estimate for " << var_name << ": " <<
                total_error << std::endl;
            }
        }

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
      libMesh::out << "Integral is " << integral << std::endl;
    }

  return 0;
}
