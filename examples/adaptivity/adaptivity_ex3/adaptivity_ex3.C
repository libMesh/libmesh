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



// <h1>Adaptivity Example 3 - Laplace Equation in the L-Shaped Domain</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example solves the Laplace equation on the classic "L-shaped"
// domain with adaptive mesh refinement.  In this case, the exact
// solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta), but
// the standard Kelly error indicator is used to estimate the error.
// The initial mesh contains three QUAD9 elements which represent the
// standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
// i.e.
// Element 0: [-1,0]x[ 0,1]
// Element 1: [ 0,1]x[ 0,1]
// Element 2: [-1,0]x[-1,0]
// The mesh is provided in the standard libMesh ASCII format file
// named "lshaped.xda".  In addition, an input file named "adaptivity_ex3.in"
// is provided which allows the user to set several parameters for
// the solution so that the problem can be re-run without a
// re-compile.  The solution technique employed is to have a
// refinement loop with a linear solve inside followed by a
// refinement of the grid and projection of the solution to the new grid
// In the final loop iteration, there is no additional
// refinement after the solve.  In the input file "adaptivity_ex3.in",
// the variable "max_r_steps" controls the number of refinement steps,
// "max_r_level" controls the maximum element refinement level, and
// "refine_percentage" / "coarsen_percentage" determine the number of
// elements which will be refined / coarsened at each step.

// LibMesh include files.
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/exact_error_estimator.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/hp_coarsentest.h"
#include "libmesh/hp_singular.h"
#include "libmesh/sibling_coupling.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/perf_log.h"
#include "libmesh/getpot.h"
#include "libmesh/exact_solution.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "libmesh/string_to_enum.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Laplace problem.  Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input.  From the
// EquationSystems object we have acess to the Mesh and
// other objects we might need.
void assemble_laplace(EquationSystems & es,
                      const std::string & system_name);


// Prototype for calculation of the exact solution.  Useful
// for setting boundary conditions.
Number exact_solution(const Point & p,
                      const Parameters &,   // EquationSystem parameters, not needed
                      const std::string &,  // sys_name, not needed
                      const std::string &); // unk_name, not needed);

// Prototype for calculation of the gradient of the exact solution.
Gradient exact_derivative(const Point & p,
                          const Parameters &,   // EquationSystems parameters, not needed
                          const std::string &,  // sys_name, not needed
                          const std::string &); // unk_name, not needed);


// These are non-const because the input file may change it,
// It is global because our exact_* functions use it.

// Set the dimensionality of the mesh
unsigned int dim = 2;

// Choose whether or not to use the singular solution
bool singularity = true;


int main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Single precision is inadequate for p refinement
  libmesh_example_requires(sizeof(Real) > 4, "--disable-singleprecision");

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Parse the input file
  GetPot input_file("adaptivity_ex3.in");

  // But allow the command line to override it.
  input_file.parse_command_line(argc, argv);

  // Read in parameters from the input file
  const unsigned int max_r_steps    = input_file("max_r_steps", 3);
  const unsigned int max_r_level    = input_file("max_r_level", 3);
  const Real refine_percentage      = input_file("refine_percentage", 0.5);
  const Real coarsen_percentage     = input_file("coarsen_percentage", 0.5);
  const unsigned int uniform_refine = input_file("uniform_refine",0);
  const std::string refine_type     = input_file("refinement_type", "h");
  const std::string approx_type     = input_file("approx_type", "LAGRANGE");
  const unsigned int approx_order   = input_file("approx_order", 1);
  const std::string element_type    = input_file("element_type", "tensor");
  const int extra_error_quadrature  = input_file("extra_error_quadrature", 0);
  const int max_linear_iterations   = input_file("max_linear_iterations", 5000);

#ifdef LIBMESH_HAVE_EXODUS_API
  const bool output_intermediate    = input_file("output_intermediate", false);
#endif

  dim = input_file("dimension", 2);
  const std::string indicator_type = input_file("indicator_type", "kelly");
  singularity = input_file("singularity", true);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // Output file for plotting the error as a function of
  // the number of degrees of freedom.
  std::string approx_name = "";
  if (element_type == "tensor")
    approx_name += "bi";
  if (approx_order == 1)
    approx_name += "linear";
  else if (approx_order == 2)
    approx_name += "quadratic";
  else if (approx_order == 3)
    approx_name += "cubic";
  else if (approx_order == 4)
    approx_name += "quartic";

  std::string output_file = approx_name;
  output_file += "_";
  output_file += refine_type;
  if (uniform_refine == 0)
    output_file += "_adaptive.m";
  else
    output_file += "_uniform.m";

  std::ofstream out (output_file.c_str());
  out << "% dofs     L2-error     H1-error" << std::endl;
  out << "e = [" << std::endl;

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system we will solve, named Laplace
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Laplace");

  // If we're doing HP refinement, then we'll need to be able to
  // evaluate data on elements' siblings in HPCoarsenTest, which means
  // we should instruct our system's DofMap to distribute that data.
  // Create and add (and keep around! the DofMap will only be holding
  // a pointer to it!) a SiblingCoupling functor to provide that
  // instruction.
  SiblingCoupling sibling_coupling;
  if (refine_type == "hp")
    system.get_dof_map().add_algebraic_ghosting_functor(sibling_coupling);

  // Read in the mesh
  if (dim == 1)
    MeshTools::Generation::build_line(mesh, 1, -1., 0.);
  else if (dim == 2)
    mesh.read("lshaped.xda");
  else
    mesh.read("lshaped3D.xda");

  // Use triangles if the config file says so
  if (element_type == "simplex")
    MeshTools::Modification::all_tri(mesh);

  // We used first order elements to describe the geometry,
  // but we may need second order elements to hold the degrees
  // of freedom
  if (approx_order > 1 || refine_type != "h")
    mesh.all_second_order();

  // Mesh Refinement object
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.refine_fraction() = refine_percentage;
  mesh_refinement.coarsen_fraction() = coarsen_percentage;
  mesh_refinement.max_h_level() = max_r_level;

  // Adds the variable "u" to "Laplace", using
  // the finite element type and order specified
  // in the config file
  system.add_variable("u", static_cast<Order>(approx_order),
                      Utility::string_to_enum<FEFamily>(approx_type));

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_laplace);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Set linear solver max iterations
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations")
    = max_linear_iterations;

  // Linear solver tolerance.
  equation_systems.parameters.set<Real>("linear solver tolerance") =
    std::pow(TOLERANCE, 2.5);

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Construct ExactSolution object and attach solution functions
  ExactSolution exact_sol(equation_systems);
  exact_sol.attach_exact_value(exact_solution);
  exact_sol.attach_exact_deriv(exact_derivative);

  // Use higher quadrature order for more accurate error results
  exact_sol.extra_quadrature_order(extra_error_quadrature);

  // A refinement loop.
  for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
    {
      libMesh::out << "Beginning Solve " << r_step << std::endl;

      // Solve the system "Laplace", just like example 2.
      system.solve();

      libMesh::out << "System has: "
                   << equation_systems.n_active_dofs()
                   << " degrees of freedom."
                   << std::endl;

      libMesh::out << "Linear solver converged at step: "
                   << system.n_linear_iterations()
                   << ", final residual: "
                   << system.final_linear_residual()
                   << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API
      // After solving the system write the solution
      // to a ExodusII-formatted plot file.
      if (output_intermediate)
        {
          std::ostringstream outfile;
          outfile << "lshaped_" << r_step << ".e";
          ExodusII_IO (mesh).write_equation_systems (outfile.str(),
                                                     equation_systems);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

      // Compute the error.
      exact_sol.compute_error("Laplace", "u");

      // Print out the error values
      libMesh::out << "L2-Error is: "
                   << exact_sol.l2_error("Laplace", "u")
                   << std::endl;
      libMesh::out << "H1-Error is: "
                   << exact_sol.h1_error("Laplace", "u")
                   << std::endl;

      // Compute any discontinuity.  There should be none.
      {
        DiscontinuityMeasure disc;
        ErrorVector disc_error;
        disc.estimate_error(system, disc_error);

        Real mean_disc_error = disc_error.mean();

        libMesh::out << "Mean discontinuity error = " << mean_disc_error << std::endl;

        // FIXME - this test fails when solving with Eigen?
#ifdef LIBMESH_ENABLE_PETSC
        libmesh_assert_less (mean_disc_error, 1e-14);
#endif
      }

      // Print to output file
      out << equation_systems.n_active_dofs() << " "
          << exact_sol.l2_error("Laplace", "u") << " "
          << exact_sol.h1_error("Laplace", "u") << std::endl;

      // Possibly refine the mesh
      if (r_step+1 != max_r_steps)
        {
          libMesh::out << "  Refining the mesh..." << std::endl;

          if (uniform_refine == 0)
            {

              // The ErrorVector is a particular StatisticsVector
              // for computing error information on a finite element mesh.
              ErrorVector error;

              if (indicator_type == "exact")
                {
                  // The ErrorEstimator class interrogates a
                  // finite element solution and assigns to each
                  // element a positive error value.
                  // This value is used for deciding which elements to
                  // refine and which to coarsen.
                  // For these simple test problems, we can use
                  // numerical quadrature of the exact error between
                  // the approximate and analytic solutions.
                  // However, for real problems, we would need an error
                  // indicator which only relies on the approximate
                  // solution.
                  ExactErrorEstimator error_estimator;

                  error_estimator.attach_exact_value(exact_solution);
                  error_estimator.attach_exact_deriv(exact_derivative);

                  // We optimize in H1 norm, the default
                  // error_estimator.error_norm = H1;

                  // Compute the error for each active element using
                  // the provided indicator.  Note in general you
                  // will need to provide an error estimator
                  // specifically designed for your application.
                  error_estimator.estimate_error (system, error);
                }
              else if (indicator_type == "patch")
                {
                  // The patch recovery estimator should give a
                  // good estimate of the solution interpolation
                  // error.
                  PatchRecoveryErrorEstimator error_estimator;

                  error_estimator.estimate_error (system, error);
                }
              else if (indicator_type == "uniform")
                {
                  // Error indication based on uniform refinement
                  // is reliable, but very expensive.
                  UniformRefinementEstimator error_estimator;

                  error_estimator.estimate_error (system, error);
                }
              else
                {
                  libmesh_assert_equal_to (indicator_type, "kelly");

                  // The Kelly error estimator is based on
                  // an error bound for the Poisson problem
                  // on linear elements, but is useful for
                  // driving adaptive refinement in many problems
                  KellyErrorEstimator error_estimator;

                  error_estimator.estimate_error (system, error);
                }

              // Write out the error distribution
              std::ostringstream ss;
              ss << r_step;
#ifdef LIBMESH_HAVE_EXODUS_API
              std::string error_output = "error_" + ss.str() + ".e";
#else
              std::string error_output = "error_" + ss.str() + ".gmv";
#endif
              error.plot_error(error_output, mesh);

              // This takes the error in error and decides which elements
              // will be coarsened or refined.  Any element within 20% of the
              // maximum error on any element will be refined, and any
              // element within 10% of the minimum error on any element might
              // be coarsened. Note that the elements flagged for refinement
              // will be refined, but those flagged for coarsening _might_ be
              // coarsened.
              mesh_refinement.flag_elements_by_error_fraction (error);

              // If we are doing adaptive p refinement, we want
              // elements flagged for that instead.
              if (refine_type == "p")
                mesh_refinement.switch_h_to_p_refinement();
              // If we are doing "matched hp" refinement, we
              // flag elements for both h and p
              else if (refine_type == "matchedhp")
                mesh_refinement.add_p_to_h_refinement();
              // If we are doing hp refinement, we
              // try switching some elements from h to p
              else if (refine_type == "hp")
                {
                  HPCoarsenTest hpselector;
                  hpselector.select_refinement(system);
                }
              // If we are doing "singular hp" refinement, we
              // try switching most elements from h to p
              else if (refine_type == "singularhp")
                {
                  // This only differs from p refinement for
                  // the singular problem
                  libmesh_assert (singularity);
                  HPSingularity hpselector;
                  // Our only singular point is at the origin
                  hpselector.singular_points.push_back(Point());
                  hpselector.select_refinement(system);
                }
              else if (refine_type != "h")
                libmesh_error_msg("Unknown refinement_type = " << refine_type);

              // This call actually refines and coarsens the flagged
              // elements.
              mesh_refinement.refine_and_coarsen_elements();
            }

          else if (uniform_refine == 1)
            {
              if (refine_type == "h" || refine_type == "hp" ||
                  refine_type == "matchedhp")
                mesh_refinement.uniformly_refine(1);
              if (refine_type == "p" || refine_type == "hp" ||
                  refine_type == "matchedhp")
                mesh_refinement.uniformly_p_refine(1);
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
  // Write out the solution
  // After solving the system write the solution
  // to a ExodusII-formatted plot file.
  ExodusII_IO (mesh).write_equation_systems ("lshaped.e",
                                             equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // Close up the output file.
  out << "];" << std::endl;
  out << "hold on" << std::endl;
  out << "plot(e(:,1), e(:,2), 'bo-');" << std::endl;
  out << "plot(e(:,1), e(:,3), 'ro-');" << std::endl;
  //    out << "set(gca,'XScale', 'Log');" << std::endl;
  //    out << "set(gca,'YScale', 'Log');" << std::endl;
  out << "xlabel('dofs');" << std::endl;
  out << "title('" << approx_name << " elements');" << std::endl;
  out << "legend('L2-error', 'H1-error');" << std::endl;
  //     out << "disp('L2-error linear fit');" << std::endl;
  //     out << "polyfit(log10(e(:,1)), log10(e(:,2)), 1)" << std::endl;
  //     out << "disp('H1-error linear fit');" << std::endl;
  //     out << "polyfit(log10(e(:,1)), log10(e(:,3)), 1)" << std::endl;
#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}




// We now define the exact solution, being careful
// to obtain an angle from atan2 in the correct
// quadrant.
Number exact_solution(const Point & p,
                      const Parameters &,  // parameters, not needed
                      const std::string &, // sys_name, not needed
                      const std::string &) // unk_name, not needed
{
  const Real x = p(0);
  const Real y = (dim > 1) ? p(1) : 0.;

  if (singularity)
    {
      // The exact solution to the singular problem,
      // u_exact = r^(2/3)*sin(2*theta/3).
      Real theta = atan2(y, x);

      // Make sure 0 <= theta <= 2*pi
      if (theta < 0)
        theta += 2. * libMesh::pi;

      // Make the 3D solution similar
      const Real z = (dim > 2) ? p(2) : 0;

      return pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
    }
  else
    {
      // The exact solution to a nonsingular problem,
      // good for testing ideal convergence rates
      const Real z = (dim > 2) ? p(2) : 0;

      return cos(x) * exp(y) * (1. - z);
    }
}





// We now define the gradient of the exact solution, again being careful
// to obtain an angle from atan2 in the correct
// quadrant.
Gradient exact_derivative(const Point & p,
                          const Parameters &,  // parameters, not needed
                          const std::string &, // sys_name, not needed
                          const std::string &) // unk_name, not needed
{
  // Gradient value to be returned.
  Gradient gradu;

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = dim > 1 ? p(1) : 0.;

  if (singularity)
    {
      // We can't compute the gradient at x=0, it is not defined.
      libmesh_assert_not_equal_to (x, 0.);

      // For convenience...
      const Real tt = 2./3.;
      const Real ot = 1./3.;

      // The value of the radius, squared
      const Real r2 = x*x + y*y;

      // The boundary value, given by the exact solution,
      // u_exact = r^(2/3)*sin(2*theta/3).
      Real theta = atan2(y, x);

      // Make sure 0 <= theta <= 2*pi
      if (theta < 0)
        theta += 2. * libMesh::pi;

      // du/dx
      gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;

      // du/dy
      if (dim > 1)
        gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;

      if (dim > 2)
        gradu(2) = 1.;
    }
  else
    {
      const Real z = (dim > 2) ? p(2) : 0;

      gradu(0) = -sin(x) * exp(y) * (1. - z);
      if (dim > 1)
        gradu(1) = cos(x) * exp(y) * (1. - z);
      if (dim > 2)
        gradu(2) = -cos(x) * exp(y);
    }

  return gradu;
}






// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_laplace(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
#ifdef LIBMESH_ENABLE_AMR
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Laplace");

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly", false);

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Laplace");

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
  UniquePtr<FEBase> fe      (FEBase::build(mesh_dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(mesh_dim, fe_type));

  // Quadrature rules for numerical integration.
  UniquePtr<QBase> qrule(fe_type.default_quadrature_rule(mesh_dim));
  UniquePtr<QBase> qface(fe_type.default_quadrature_rule(mesh_dim-1));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (qrule.get());
  fe_face->attach_quadrature_rule (qface.get());

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties or forcing functions at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  // For this simple problem we usually only need them on element
  // boundaries.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

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
  // example 3 for a discussion of the element iterators.  Here we use
  // the const_active_local_elem_iterator to indicate we only want
  // to loop over elements that are assigned to the local processor
  // which are "active" in the sense of AMR.  This allows each
  // processor to compute its components of the global matrix for
  // active elements while ignoring parent elements which have been
  // refined.
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
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Stop logging the shape function initialization.
      // If you forget to stop logging an event the PerfLog
      // object will probably catch the error and abort.
      perf_log.pop("elem init");

      // Now we will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      //
      // Now start logging the element matrix computation
      perf_log.push ("Ke");

      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        for (std::size_t i=0; i<dphi.size(); i++)
          for (std::size_t j=0; j<dphi.size(); j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      // We need a forcing function to make the 1D case interesting
      if (mesh_dim == 1)
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          {
            Real x = q_point[qp](0);
            Real f = singularity ? sqrt(3.)/9.*pow(-x, -4./3.) :
              cos(x);
            for (std::size_t i=0; i<dphi.size(); ++i)
              Fe(i) += JxW[qp]*phi[i][qp]*f;
          }

      // Stop logging the matrix computation
      perf_log.pop ("Ke");


      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.
      //
      // This approach adds the L2 projection of the boundary
      // data in penalty form to the weak statement.  This is
      // a more generic approach for applying Dirichlet BCs
      // which is applicable to non-Lagrange finite element
      // discretizations.
      {
        // Start logging the boundary condition computation
        perf_log.push ("BCs");

        // The penalty value.
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              fe_face->reinit(elem, s);

              for (unsigned int qp=0; qp<qface->n_points(); qp++)
                {
                  const Number value = exact_solution (qface_points[qp],
                                                       es.parameters,
                                                       "null",
                                                       "void");

                  // RHS contribution
                  for (std::size_t i=0; i<psi.size(); i++)
                    Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];

                  // Matrix contribution
                  for (std::size_t i=0; i<psi.size(); i++)
                    for (std::size_t j=0; j<psi.size(); j++)
                      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                }
            }

        // Stop logging the boundary condition computation
        perf_log.pop ("BCs");
      }


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

      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.pop ("matrix insertion");
    }

  // That's it.  We don't need to do anything else to the
  // PerfLog.  When it goes out of scope (at this function return)
  // it will print its log to the screen. Pretty easy, huh?
#endif // #ifdef LIBMESH_ENABLE_AMR
}
