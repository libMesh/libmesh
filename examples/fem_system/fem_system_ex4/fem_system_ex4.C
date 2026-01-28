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



// <h1>FEMSystem Example 4 - Mixed-dimension heat transfer equation
// with FEMSystem</h1>
// \author Roy Stogner
// \date 2015
//
// This example shows how elements of multiple dimensions can be
// linked and computed upon simultaneously within the
// DifferentiableSystem class framework

// C++ includes
#include <iomanip>
#include <memory>

// Basic include files
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/exact_solution.h"
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/parsed_function.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_norm_type.h"

// The systems and solvers we may use
#include "heatsystem.h"
#include "libmesh/diff_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // This doesn't converge without at least double precision
  libmesh_example_requires(sizeof(Real) > 4, "--disable-singleprecision");

  // Parse the input file
  GetPot infile("fem_system_ex4.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the input file
  const Real global_tolerance          = infile("global_tolerance", 0.);
  const unsigned int nelem_target      = infile("n_elements", 400);
  const Real deltat                    = infile("deltat", 0.005);
  const unsigned int coarsegridsize    = infile("coarsegridsize", 20);
  const unsigned int coarserefinements = infile("coarserefinements", 0);
  const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 10);
  const unsigned int dim               = infile("dimension", 2);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // We have only defined 2 and 3 dimensional problems
  libmesh_assert (dim == 2 || dim == 3);

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // And an object to refine it
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.coarsen_by_parents() = true;
  mesh_refinement.absolute_global_tolerance() = global_tolerance;
  mesh_refinement.nelem_target() = nelem_target;
  mesh_refinement.refine_fraction() = 0.3;
  mesh_refinement.coarsen_fraction() = 0.3;
  mesh_refinement.coarsen_threshold() = 0.1;

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square or cube.  We crop the domain at y=2/3 to allow
  // for a homogeneous Neumann BC in our benchmark there.
  boundary_id_type bcid = 3; // +y in 3D
  if (dim == 2)
    {
      MeshTools::Generation::build_square
        (mesh,
         coarsegridsize,
         coarsegridsize*2/3, // int arithmetic best we can do here
         0., 1.,
         0., 2./3.,
         QUAD9);
      bcid = 2; // +y in 2D
    }
  else if (dim == 3)
    {
      MeshTools::Generation::build_cube
        (mesh,
         coarsegridsize,
         coarsegridsize*2/3,
         coarsegridsize,
         0., 1.,
         0., 2./3.,
         0., 1.,
         HEX27);
    }

  {
    // Add boundary elements corresponding to the +y boundary of our
    // volume mesh
    std::set<boundary_id_type> bcids;
    bcids.insert(bcid);
    mesh.get_boundary_info().add_elements(bcids, mesh);
    mesh.prepare_for_use();
  }

  // To work around ExodusII file format limitations, we need elements
  // of different dimensionality to belong to different subdomains.
  // Our interior elements defaulted to subdomain id 0, so we'll set
  // boundary elements to subdomain 1.
  for (auto & elem : mesh.element_ptr_range())
    if (elem->dim() < dim)
      elem->subdomain_id() = 1;

  // Make sure the mesh knows we added new subdomains.
  mesh.cache_elem_data();

  mesh_refinement.uniformly_refine(coarserefinements);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system "Heat" and its variables.
  HeatSystem & system =
    equation_systems.add_system<HeatSystem> ("Heat");

  // Solve this as a steady system
  system.time_solver = std::make_unique<SteadySolver>(system);

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver & solver = *(system.time_solver->diff_solver().get());
  solver.quiet = infile("solver_quiet", true);
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations = infile("max_nonlinear_iterations", 15);
  solver.relative_step_tolerance = infile("relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance = infile("relative_residual_tolerance", 0.0);
  solver.absolute_residual_tolerance = infile("absolute_residual_tolerance", 0.0);

  // And the linear solver options
  solver.max_linear_iterations = infile("max_linear_iterations", 50000);
  solver.initial_linear_tolerance = infile("initial_linear_tolerance", 1.e-3);

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Adaptively solve the steady solution
  unsigned int a_step = 0;
  for (; a_step != max_adaptivesteps; ++a_step)
    {
      system.solve();

      system.postprocess();

      ErrorVector error;

      std::unique_ptr<ErrorEstimator> error_estimator;

      // To solve to a tolerance in this problem we
      // need a better estimator than Kelly
      if (global_tolerance != 0.)
        {
          // We can't adapt to both a tolerance and a mesh
          // size at once
          libmesh_assert_equal_to (nelem_target, 0);

          auto u = std::make_unique<UniformRefinementEstimator>();

          // The lid-driven cavity problem isn't in H1, so
          // lets estimate L2 error
          u->error_norm = L2;
          error_estimator = std::move(u);
        }
      else
        {
          // If we aren't adapting to a tolerance we need a
          // target mesh size
          libmesh_assert_greater (nelem_target, 0);

          // Kelly is a lousy estimator to use for a problem
          // not in H1 - if we were doing more than a few
          // timesteps we'd need to turn off or limit the
          // maximum level of our adaptivity eventually
          error_estimator = std::make_unique<KellyErrorEstimator>();
        }

      error_estimator->estimate_error(system, error);

      // Print out status at each adaptive step.
      Real global_error = error.l2_norm();
      libMesh::out << "Adaptive step "
                   << a_step
                   << ": "
                   << std::endl;

      if (global_tolerance != 0.)
        libMesh::out << "Global_error = "
                     << global_error
                     << std::endl;

      if (global_tolerance != 0.)
        libMesh::out << "Worst element error = "
                     << error.maximum()
                     << ", mean = "
                     << error.mean()
                     << std::endl;

      if (global_tolerance != 0.)
        {
          // If we've reached our desired tolerance, we
          // don't need any more adaptive steps
          if (global_error < global_tolerance)
            break;
          mesh_refinement.flag_elements_by_error_tolerance(error);
        }
      else
        {
          // If flag_elements_by_nelem_target returns true, this
          // should be our last adaptive step.
          if (mesh_refinement.flag_elements_by_nelem_target(error))
            {
              mesh_refinement.refine_and_coarsen_elements();
              equation_systems.reinit();
              a_step = max_adaptivesteps;
              break;
            }
        }

      // Carry out the adaptive mesh refinement/coarsening
      mesh_refinement.refine_and_coarsen_elements();
      equation_systems.reinit();

      libMesh::out << "Refined mesh to "
                   << mesh.n_active_elem()
                   << " active elements and "
                   << equation_systems.n_active_dofs()
                   << " active dofs."
                   << std::endl;
    }
  // Do one last solve if necessary
  if (a_step == max_adaptivesteps)
    {
      system.solve();

      system.postprocess();
    }


#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems
    ("out.e", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

#ifdef LIBMESH_HAVE_GMV
  GMVIO(mesh).write_equation_systems
    ("out.gmv", equation_systems);
#endif // #ifdef LIBMESH_HAVE_GMV

#ifdef LIBMESH_HAVE_FPARSER
  // Check that we got close to the analytic solution
  ExactSolution exact_sol(equation_systems);
  const std::string exact_str = (dim == 2) ?
    "sin(pi*x)*sin(pi*y)" : "sin(pi*x)*sin(pi*y)*sin(pi*z)";
  ParsedFunction<Number> exact_func(exact_str);
  exact_sol.attach_exact_value(0, &exact_func);
  exact_sol.compute_error("Heat", "T");

  Real err = exact_sol.l2_error("Heat", "T");

  // Print out the error value
  libMesh::out << "L2-Error is: " << err << std::endl;

  libmesh_assert_less(err, 2e-3);

#endif // #ifdef LIBMESH_HAVE_FPARSER

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
