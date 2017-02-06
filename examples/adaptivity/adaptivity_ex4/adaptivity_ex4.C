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

// <h1>Adaptivity Example 4 - Biharmonic Equation</h1>
// \author Benjamin S. Kirk
// \date 2004
//
// This example solves the Biharmonic equation on a square or cube,
// using a Galerkin formulation with C1 elements approximating the
// H^2_0 function space.
// The initial mesh contains two TRI6, one QUAD9 or one HEX27
// An input file named "ex15.in"
// is provided which allows the user to set several parameters for
// the solution so that the problem can be re-run without a
// re-compile.  The solution technique employed is to have a
// refinement loop with a linear solve inside followed by a
// refinement of the grid and projection of the solution to the new grid
// In the final loop iteration, there is no additional
// refinement after the solve.  In the input file "ex15.in", the variable
// "max_r_steps" controls the number of refinement steps, and
// "max_r_level" controls the maximum element refinement level.

// LibMesh include files.
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/fourth_error_estimators.h"
#include "libmesh/getpot.h"
#include "libmesh/exact_solution.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"
#include "libmesh/perf_log.h"
#include "libmesh/string_to_enum.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Biharmonic problem.  Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input.  From the
// EquationSystems object we have acess to the Mesh and
// other objects we might need.
void assemble_biharmonic(EquationSystems & es,
                         const std::string & system_name);


// Prototypes for calculation of the exact solution.  Necessary
// for setting boundary conditions.
Number exact_1D_solution(const Point & p,
                         const Parameters &,
                         const std::string &,
                         const std::string &);

Number exact_2D_solution(const Point & p,
                         const Parameters &,   // parameters, not needed
                         const std::string &,  // sys_name, not needed
                         const std::string &); // unk_name, not needed);

Number exact_3D_solution(const Point & p,
                         const Parameters &,
                         const std::string &,
                         const std::string &);

// Prototypes for calculation of the gradient of the exact solution.
// Necessary for setting boundary conditions in H^2_0 and testing
// H^1 convergence of the solution
Gradient exact_1D_derivative(const Point & p,
                             const Parameters &,
                             const std::string &,
                             const std::string &);

Gradient exact_2D_derivative(const Point & p,
                             const Parameters &,
                             const std::string &,
                             const std::string &);

Gradient exact_3D_derivative(const Point & p,
                             const Parameters &,
                             const std::string &,
                             const std::string &);

Tensor exact_1D_hessian(const Point & p,
                        const Parameters &,
                        const std::string &,
                        const std::string &);

Tensor exact_2D_hessian(const Point & p,
                        const Parameters &,
                        const std::string &,
                        const std::string &);

Tensor exact_3D_hessian(const Point & p,
                        const Parameters &,
                        const std::string &,
                        const std::string &);

Number forcing_function_1D(const Point & p);

Number forcing_function_2D(const Point & p);

Number forcing_function_3D(const Point & p);

// Pointers to dimension-independent functions
Number (*exact_solution)(const Point & p,
                         const Parameters &,
                         const std::string &,
                         const std::string &);

Gradient (*exact_derivative)(const Point & p,
                             const Parameters &,
                             const std::string &,
                             const std::string &);

Tensor (*exact_hessian)(const Point & p,
                        const Parameters &,
                        const std::string &,
                        const std::string &);

Number (*forcing_function)(const Point & p);



int main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Adaptive constraint calculations for fine Hermite elements seems
  // to require half-decent precision
#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  libmesh_example_requires(false, "double precision");
#endif

  // This example requires Adaptive Mesh Refinement support
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // This example requires second derivative calculation support
#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
  libmesh_example_requires(false, "--enable-second");
#else

  // Parse the input file
  GetPot input_file("adaptivity_ex4.in");

  // But allow the command line to override it.
  input_file.parse_command_line(argc, argv);

  // Read in parameters from the input file
  const unsigned int max_r_level = input_file("max_r_level", 10);
  const unsigned int max_r_steps = input_file("max_r_steps", 4);
  const std::string approx_type  = input_file("approx_type", "HERMITE");
  const std::string approx_order_string = input_file("approx_order", "THIRD");
  const unsigned int uniform_refine = input_file("uniform_refine", 0);
  const Real refine_percentage = input_file("refine_percentage", 0.5);
  const Real coarsen_percentage = input_file("coarsen_percentage", 0.5);
  const unsigned int dim = input_file("dimension", 2);
  const unsigned int max_linear_iterations = input_file("max_linear_iterations", 10000);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // Currently only the Hermite cubics give a 1D or 3D C^1 basis
  libmesh_assert (dim == 2 || approx_type == "HERMITE");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Output file for plotting the error
  std::string output_file = "";

  if (dim == 1)
    output_file += "1D_";
  else if (dim == 2)
    output_file += "2D_";
  else if (dim == 3)
    output_file += "3D_";

  if (approx_type == "HERMITE")
    output_file += "hermite_";
  else if (approx_type == "SECOND")
    output_file += "reducedclough_";
  else
    output_file += "clough_";

  if (uniform_refine == 0)
    output_file += "adaptive";
  else
    output_file += "uniform";

#ifdef LIBMESH_HAVE_EXODUS_API
  // If we have Exodus, use the same base output filename
  std::string exd_file = output_file;
  exd_file += ".e";
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  output_file += ".m";

  std::ofstream out (output_file.c_str());
  out << "% dofs     L2-error     H1-error      H2-error\n"
      << "e = [\n";

  // Set up the dimension-dependent coarse mesh and solution
  if (dim == 1)
    {
      MeshTools::Generation::build_line(mesh, 2);
      exact_solution = &exact_1D_solution;
      exact_derivative = &exact_1D_derivative;
      exact_hessian = &exact_1D_hessian;
      forcing_function = &forcing_function_1D;
    }

  if (dim == 2)
    {
      MeshTools::Generation::build_square(mesh, 2, 2);
      exact_solution = &exact_2D_solution;
      exact_derivative = &exact_2D_derivative;
      exact_hessian = &exact_2D_hessian;
      forcing_function = &forcing_function_2D;
    }
  else if (dim == 3)
    {
      MeshTools::Generation::build_cube(mesh, 2, 2, 2);
      exact_solution = &exact_3D_solution;
      exact_derivative = &exact_3D_derivative;
      exact_hessian = &exact_3D_hessian;
      forcing_function = &forcing_function_3D;
    }

  // Convert the mesh to second order: necessary for computing with
  // Clough-Tocher elements, useful for getting slightly less
  // broken visualization output with Hermite elements
  mesh.all_second_order();

  // Convert it to triangles if necessary
  if (approx_type != "HERMITE")
    MeshTools::Modification::all_tri(mesh);

  // Mesh Refinement object
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.refine_fraction() = refine_percentage;
  mesh_refinement.coarsen_fraction() = coarsen_percentage;
  mesh_refinement.max_h_level() = max_r_level;

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Biharmonic"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Biharmonic");

  Order approx_order = approx_type == "SECOND" ? SECOND :
    Utility::string_to_enum<Order>(approx_order_string);

  // Adds the variable "u" to "Biharmonic".  "u" will be approximated
  // using Hermite tensor product squares or Clough-Tocher triangles

  if (approx_type == "HERMITE")
    system.add_variable("u", approx_order, HERMITE);
  else if (approx_type == "SECOND")
    system.add_variable("u", SECOND, CLOUGH);
  else if (approx_type == "CLOUGH")
    system.add_variable("u", approx_order, CLOUGH);
  else
    libmesh_error_msg("Invalid approx_type = " << approx_type);

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_biharmonic);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Set linear solver max iterations
  equation_systems.parameters.set<unsigned int>
    ("linear solver maximum iterations") = max_linear_iterations;

  // Linear solver tolerance.
  equation_systems.parameters.set<Real>
    ("linear solver tolerance") = TOLERANCE * TOLERANCE;

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Construct ExactSolution object and attach function to compute exact solution
  ExactSolution exact_sol(equation_systems);
  exact_sol.attach_exact_value(exact_solution);
  exact_sol.attach_exact_deriv(exact_derivative);
  exact_sol.attach_exact_hessian(exact_hessian);

  // Construct zero solution object, useful for computing solution norms
  // Attaching "zero_solution" functions is unnecessary
  ExactSolution zero_sol(equation_systems);

  // A refinement loop.
  for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
    {
      mesh.print_info();
      equation_systems.print_info();

      libMesh::out << "Beginning Solve " << r_step << std::endl;

      // Solve the system "Biharmonic", just like example 2.
      system.solve();

      libMesh::out << "Linear solver converged at step: "
                   << system.n_linear_iterations()
                   << ", final residual: "
                   << system.final_linear_residual()
                   << std::endl;

      // Compute the error.
      exact_sol.compute_error("Biharmonic", "u");
      // Compute the norm.
      zero_sol.compute_error("Biharmonic", "u");

      // Print out the error values
      libMesh::out << "L2-Norm is: "
                   << zero_sol.l2_error("Biharmonic", "u")
                   << std::endl;
      libMesh::out << "H1-Norm is: "
                   << zero_sol.h1_error("Biharmonic", "u")
                   << std::endl;
      libMesh::out << "H2-Norm is: "
                   << zero_sol.h2_error("Biharmonic", "u")
                   << std::endl
                   << std::endl;
      libMesh::out << "L2-Error is: "
                   << exact_sol.l2_error("Biharmonic", "u")
                   << std::endl;
      libMesh::out << "H1-Error is: "
                   << exact_sol.h1_error("Biharmonic", "u")
                   << std::endl;
      libMesh::out << "H2-Error is: "
                   << exact_sol.h2_error("Biharmonic", "u")
                   << std::endl
                   << std::endl;

      // Print to output file
      out << equation_systems.n_active_dofs() << " "
          << exact_sol.l2_error("Biharmonic", "u") << " "
          << exact_sol.h1_error("Biharmonic", "u") << " "
          << exact_sol.h2_error("Biharmonic", "u") << std::endl;

      // Possibly refine the mesh
      if (r_step+1 != max_r_steps)
        {
          libMesh::out << "  Refining the mesh..." << std::endl;

          if (uniform_refine == 0)
            {
              ErrorVector error;
              LaplacianErrorEstimator error_estimator;

              error_estimator.estimate_error(system, error);
              mesh_refinement.flag_elements_by_elem_fraction (error);

              libMesh::out << "Mean Error: " << error.mean() << std::endl;
              libMesh::out << "Error Variance: " << error.variance() << std::endl;

              mesh_refinement.refine_and_coarsen_elements();
            }
          else
            {
              mesh_refinement.uniformly_refine(1);
            }

          // This call reinitializes the EquationSystems object for
          // the newly refined mesh.  One of the steps in the
          // reinitialization is projecting the solution,
          // old_solution, etc... vectors from the old mesh to
          // the current one.
          equation_systems.reinit ();
        }
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  // After solving the system write the solution
  // to a ExodusII-formatted plot file.
  ExodusII_IO (mesh).write_equation_systems (exd_file,
                                             equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // Close up the output file.
  out << "];\n"
      << "hold on\n"
      << "loglog(e(:,1), e(:,2), 'bo-');\n"
      << "loglog(e(:,1), e(:,3), 'ro-');\n"
      << "loglog(e(:,1), e(:,4), 'go-');\n"
      << "xlabel('log(dofs)');\n"
      << "ylabel('log(error)');\n"
      << "title('C1 " << approx_type << " elements');\n"
      << "legend('L2-error', 'H1-error', 'H2-error');\n";

  // All done.
  return 0;
#endif // #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
#endif // #ifndef LIBMESH_ENABLE_AMR
}



Number exact_1D_solution(const Point & p,
                         const Parameters &,  // parameters, not needed
                         const std::string &, // sys_name, not needed
                         const std::string &) // unk_name, not needed
{
  // x coordinate in space
  const Real x = p(0);

  // analytic solution value
  return 256.*(x-x*x)*(x-x*x);
}


// We now define the gradient of the exact solution
Gradient exact_1D_derivative(const Point & p,
                             const Parameters &,  // parameters, not needed
                             const std::string &, // sys_name, not needed
                             const std::string &) // unk_name, not needed
{
  // x coordinate in space
  const Real x = p(0);

  // First derivatives to be returned.
  Gradient gradu;

  gradu(0) = 256.*2.*(x-x*x)*(1-2*x);

  return gradu;
}


// We now define the hessian of the exact solution
Tensor exact_1D_hessian(const Point & p,
                        const Parameters &,  // parameters, not needed
                        const std::string &, // sys_name, not needed
                        const std::string &) // unk_name, not needed
{
  // Second derivatives to be returned.
  Tensor hessu;

  // x coordinate in space
  const Real x = p(0);

  hessu(0,0) = 256.*2.*(1-6.*x+6.*x*x);

  return hessu;
}



Number forcing_function_1D(const Point &)
{
  // Equals laplacian(laplacian(u)), u'''' in 1D
  return 256. * 2. * 12.;
}


Number exact_2D_solution(const Point & p,
                         const Parameters &,  // parameters, not needed
                         const std::string &, // sys_name, not needed
                         const std::string &) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);

  // analytic solution value
  return 256.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y);
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative(const Point & p,
                             const Parameters &,  // parameters, not needed
                             const std::string &, // sys_name, not needed
                             const std::string &) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);

  // First derivatives to be returned.
  Gradient gradu;

  gradu(0) = 256.*2.*(x-x*x)*(1-2*x)*(y-y*y)*(y-y*y);
  gradu(1) = 256.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1-2*y);

  return gradu;
}


// We now define the hessian of the exact solution
Tensor exact_2D_hessian(const Point & p,
                        const Parameters &,  // parameters, not needed
                        const std::string &, // sys_name, not needed
                        const std::string &) // unk_name, not needed
{
  // Second derivatives to be returned.
  Tensor hessu;

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);

  hessu(0,0) = 256.*2.*(1-6.*x+6.*x*x)*(y-y*y)*(y-y*y);
  hessu(0,1) = 256.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y);
  hessu(1,1) = 256.*2.*(x-x*x)*(x-x*x)*(1.-6.*y+6.*y*y);

  // Hessians are always symmetric
  hessu(1,0) = hessu(0,1);
  return hessu;
}



Number forcing_function_2D(const Point & p)
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);

  // Equals laplacian(laplacian(u))
  return 256. * 8. * (3.*((y-y*y)*(y-y*y)+(x-x*x)*(x-x*x))
                      + (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y));
}



Number exact_3D_solution(const Point & p,
                         const Parameters &,  // parameters, not needed
                         const std::string &, // sys_name, not needed
                         const std::string &) // unk_name, not needed
{
  // xyz coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  // analytic solution value
  return 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
}


Gradient exact_3D_derivative(const Point & p,
                             const Parameters &,  // parameters, not needed
                             const std::string &, // sys_name, not needed
                             const std::string &) // unk_name, not needed
{
  // First derivatives to be returned.
  Gradient gradu;

  // xyz coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  gradu(0) = 4096.*2.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
  gradu(1) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
  gradu(2) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);

  return gradu;
}


// We now define the hessian of the exact solution
Tensor exact_3D_hessian(const Point & p,
                        const Parameters &,  // parameters, not needed
                        const std::string &, // sys_name, not needed
                        const std::string &) // unk_name, not needed
{
  // Second derivatives to be returned.
  Tensor hessu;

  // xyz coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  hessu(0,0) = 4096.*(2.-12.*x+12.*x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
  hessu(0,1) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
  hessu(0,2) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);
  hessu(1,1) = 4096.*(x-x*x)*(x-x*x)*(2.-12.*y+12.*y*y)*(z-z*z)*(z-z*z);
  hessu(1,2) = 4096.*4.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(1.-2.*z);
  hessu(2,2) = 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(2.-12.*z+12.*z*z);

  // Hessians are always symmetric
  hessu(1,0) = hessu(0,1);
  hessu(2,0) = hessu(0,2);
  hessu(2,1) = hessu(1,2);

  return hessu;
}



Number forcing_function_3D(const Point & p)
{
  // xyz coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  // Equals laplacian(laplacian(u))
  return 4096. * 8. * (3.*((y-y*y)*(y-y*y)*(x-x*x)*(x-x*x) +
                           (z-z*z)*(z-z*z)*(x-x*x)*(x-x*x) +
                           (z-z*z)*(z-z*z)*(y-y*y)*(y-y*y)) +
                       (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y)*(z-z*z)*(z-z*z) +
                       (1.-6.*x+6.*x*x)*(1.-6.*z+6.*z*z)*(y-y*y)*(y-y*y) +
                       (1.-6.*y+6.*y*y)*(1.-6.*z+6.*z*z)*(x-x*x)*(x-x*x));
}



// We now define the matrix assembly function for the
// Biharmonic system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_biharmonic(EquationSystems & es,
                         const std::string & libmesh_dbg_var(system_name))
{
#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Biharmonic");

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly", false);

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Biharmonic");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  // Quadrature rule for numerical integration.
  // With 2D triangles, the Clough quadrature rule puts a Gaussian
  // quadrature rule on each of the 3 subelements
  UniquePtr<QBase> qrule(fe_type.default_quadrature_rule(dim));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (qrule.get());

  // Declare a special finite element object for
  // boundary integration.
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires another quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  // In 1D, the Clough and Gauss quadrature rules are identical.
  UniquePtr<QBase> qface(fe_type.default_quadrature_rule(dim-1));

  // Tell the finte element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (qface.get());

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function second derivatives evaluated at the
  // quadrature points.  Note that for the simple biharmonic, shape
  // function first derivatives are unnecessary.
  const std::vector<std::vector<RealTensor> > & d2phi = fe->get_d2phi();

  // For efficiency we will compute shape function laplacians n times,
  // not n^2
  std::vector<Real> shape_laplacian;

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe". More detail is in example 3.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.  See
  // example 3 for a discussion of the element iterators.

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Start logging the shape function initialization.
      // This is done through a simple function call with
      // the name of the event to log.
      perf_log.push("elem init");

      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Make sure there is enough room in this cache
      shape_laplacian.resize(dof_indices.size());

      // Stop logging the shape function initialization.
      // If you forget to stop logging an event the PerfLog
      // object will probably catch the error and abort.
      perf_log.pop("elem init");

      // Now we will build the element matrix.  This involves
      // a double loop to integrate laplacians of the test funcions
      // (i) against laplacians of the trial functions (j).
      //
      // This step is why we need the Clough-Tocher elements -
      // these C1 differentiable elements have square-integrable
      // second derivatives.
      //
      // Now start logging the element matrix computation
      perf_log.push ("Ke");

      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          for (std::size_t i=0; i<phi.size(); i++)
            {
              shape_laplacian[i] = d2phi[i][qp](0,0);
              if (dim > 1)
                shape_laplacian[i] += d2phi[i][qp](1,1);
              if (dim == 3)
                shape_laplacian[i] += d2phi[i][qp](2,2);
            }
          for (std::size_t i=0; i<phi.size(); i++)
            for (std::size_t j=0; j<phi.size(); j++)
              Ke(i,j) += JxW[qp]*
                shape_laplacian[i]*shape_laplacian[j];
        }

      // Stop logging the matrix computation
      perf_log.pop ("Ke");


      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.  Note that this is a fourth-order
      // problem: Dirichlet boundary conditions include *both*
      // boundary values and boundary normal fluxes.
      {
        // Start logging the boundary condition computation
        perf_log.push ("BCs");

        // The penalty values, for solution boundary trace and flux.
        const Real penalty = 1e10;
        const Real penalty2 = 1e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector<Real> > & phi_face =
                fe_face->get_phi();

              // The value of the shape function derivatives at the
              // quadrature points.
              const std::vector<std::vector<RealGradient> > & dphi_face =
                fe_face->get_dphi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point> & qface_point = fe_face->get_xyz();
              const std::vector<Point> & face_normals = fe_face->get_normals();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, s);

              // Loop over the face quagrature points for integration.
              for (unsigned int qp=0; qp<qface->n_points(); qp++)
                {
                  // The boundary value.
                  Number value = exact_solution(qface_point[qp],
                                                es.parameters, "null",
                                                "void");
                  Gradient flux = exact_2D_derivative(qface_point[qp],
                                                      es.parameters,
                                                      "null", "void");

                  // Matrix contribution of the L2 projection.
                  // Note that the basis function values are
                  // integrated against test function values while
                  // basis fluxes are integrated against test function
                  // fluxes.
                  for (std::size_t i=0; i<phi_face.size(); i++)
                    for (std::size_t j=0; j<phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp] *
                        (penalty * phi_face[i][qp] *
                         phi_face[j][qp] + penalty2
                         * (dphi_face[i][qp] *
                            face_normals[qp]) *
                         (dphi_face[j][qp] *
                          face_normals[qp]));

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (std::size_t i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp] *
                      (penalty * value * phi_face[i][qp]
                       + penalty2 *
                       (flux * face_normals[qp])
                       * (dphi_face[i][qp]
                          * face_normals[qp]));

                }
            }

        // Stop logging the boundary condition computation
        perf_log.pop ("BCs");
      }

      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        for (std::size_t i=0; i<phi.size(); i++)
          Fe(i) += JxW[qp]*phi[i][qp]*forcing_function(q_point[qp]);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.push ("matrix insertion");

      dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

      // Stop logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.pop ("matrix insertion");
    }

  // That's it.  We don't need to do anything else to the
  // PerfLog.  When it goes out of scope (at this function return)
  // it will print its log to the screen. Pretty easy, huh?

#else

#endif // #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#endif // #ifdef LIBMESH_ENABLE_AMR
}
