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



// <h1>Introduction Example 4 - Solving a 1D, 2D or 3D Poisson Problem in Parallel</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the fourth example program.  It builds on
// the third example program by showing how to formulate
// the code in a dimension-independent way.  Very minor
// changes to the example will allow the problem to be
// solved in one, two or three dimensions.
//
// This example will also introduce the PerfLog class
// as a way to monitor your code's performance.  We will
// use it to instrument the matrix assembly code and look
// for bottlenecks where we should focus optimization efforts.
//
// This example also shows how to extend example 3 to run in
// parallel.  Notice how little has changed!


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input.  From the
// EquationSystems object we have acess to the Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems & es,
                      const std::string & system_name);

// Exact solution function prototype.
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.);

// Define a wrapper for exact_solution that will be needed below
void exact_solution_wrapper (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
  output(0) = exact_solution(p(0),
                             (LIBMESH_DIM>1)?p(1):0,
                             (LIBMESH_DIM>2)?p(2):0);
}

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);

  // Declare a performance log for the main program
  // PerfLog perf_main("Main Program");

  // Create a GetPot object to parse the command line
  GetPot command_line (argc, argv);

  // Check for proper calling arguments.
  if (argc < 3)
    {
      // This handy function will print the file name, line number,
      // specified message, and then throw an exception.
      libmesh_error_msg("Usage:\n" << "\t " << argv[0] << " -d 2(3)" << " -n 15");
    }

  // Brief message to the user regarding the program name
  // and command line arguments.
  else
    {
      libMesh::out << "Running " << argv[0];

      for (int i=1; i<argc; i++)
        libMesh::out << " " << argv[i];

      libMesh::out << std::endl << std::endl;
    }

  // Read problem dimension from command line.  Use int
  // instead of unsigned since the GetPot overload is ambiguous
  // otherwise.
  int dim = 2;
  if (command_line.search(1, "-d"))
    dim = command_line.next(dim);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // Create a mesh with user-defined dimension.
  // Read number of elements from command line
  int ps = 15;
  if (command_line.search(1, "-n"))
    ps = command_line.next(ps);

  // Read FE order from command line
  std::string order = "SECOND";
  if (command_line.search(2, "-Order", "-o"))
    order = command_line.next(order);

  // Read FE Family from command line
  std::string family = "LAGRANGE";
  if (command_line.search(2, "-FEFamily", "-f"))
    family = command_line.next(family);

  // Cannot use discontinuous basis.
  if ((family == "MONOMIAL") || (family == "XYZ"))
    libmesh_error_msg("ex4 currently requires a C^0 (or higher) FE basis.");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements in 2D, or Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.

  Real halfwidth = dim > 1 ? 1. : 0.;
  Real halfheight = dim > 2 ? 1. : 0.;

  if ((family == "LAGRANGE") && (order == "FIRST"))
    {
      // No reason to use high-order geometric elements if we are
      // solving with low-order finite elements.
      MeshTools::Generation::build_cube (mesh,
                                         ps,
                                         (dim>1) ? ps : 0,
                                         (dim>2) ? ps : 0,
                                         -1., 1.,
                                         -halfwidth, halfwidth,
                                         -halfheight, halfheight,
                                         (dim==1)    ? EDGE2 :
                                         ((dim == 2) ? QUAD4 : HEX8));
    }

  else
    {
      MeshTools::Generation::build_cube (mesh,
                                         ps,
                                         (dim>1) ? ps : 0,
                                         (dim>2) ? ps : 0,
                                         -1., 1.,
                                         -halfwidth, halfwidth,
                                         -halfheight, halfheight,
                                         (dim==1)    ? EDGE3 :
                                         ((dim == 2) ? QUAD9 : HEX27));
    }


  // Print information about the mesh to the screen.
  mesh.print_info();


  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Poisson"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Poisson");


  // Add the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  unsigned int u_var = system.add_variable("u",
                                           Utility::string_to_enum<Order>   (order),
                                           Utility::string_to_enum<FEFamily>(family));

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_poisson);

  // Construct a Dirichlet boundary condition object

  // Indicate which boundary IDs we impose the BC on
  // We either build a line, a square or a cube, and
  // here we indicate the boundaries IDs in each case
  std::set<boundary_id_type> boundary_ids;
  // the dim==1 mesh has two boundaries with IDs 0 and 1
  boundary_ids.insert(0);
  boundary_ids.insert(1);
  // the dim==2 mesh has four boundaries with IDs 0, 1, 2 and 3
  if (dim>=2)
    {
      boundary_ids.insert(2);
      boundary_ids.insert(3);
    }
  // the dim==3 mesh has four boundaries with IDs 0, 1, 2, 3, 4 and 5
  if (dim==3)
    {
      boundary_ids.insert(4);
      boundary_ids.insert(5);
    }

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables(1);
  variables[0] = u_var;

  // Create an AnalyticFunction object that we use to project the BC
  // This function just calls the function exact_solution via exact_solution_wrapper
  AnalyticFunction<> exact_solution_object(exact_solution_wrapper);

  // In general, when reusing a system-indexed exact solution, we want
  // to use the default system-ordering constructor for
  // DirichletBoundary, so we demonstrate that here.  In this case,
  // though, we have only one variable, so system- and local-
  // orderings are the same.
  DirichletBoundary dirichlet_bc
    (boundary_ids, variables, exact_solution_object);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);


  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();
  mesh.print_info();

  // Solve the system "Poisson", just like example 2.
  system.solve();

  // After solving the system write the solution
  // to a GMV-formatted plot file.
  if (dim == 1)
    {
      GnuPlotIO plot(mesh, "Introduction Example 4, 1D", GnuPlotIO::GRID_ON);
      plot.write_equation_systems("gnuplot_script", equation_systems);
    }
#ifdef LIBMESH_HAVE_EXODUS_API
  else
    {
      ExodusII_IO (mesh).write_equation_systems ((dim == 3) ?
                                                 "out_3.e" : "out_2.e", equation_systems);
    }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}




// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions.
void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Poisson");

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Poisson");

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

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);

  // Tell the finte element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

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

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

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

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.  See example 3 for a discussion of the
  // element iterators.
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
      // We have split the numeric integration into two loops
      // so that we can log the matrix and right-hand-side
      // computation seperately.
      //
      // Now start logging the element matrix computation
      perf_log.push ("Ke");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (std::size_t i=0; i<phi.size(); i++)
          for (std::size_t j=0; j<phi.size(); j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);


      // Stop logging the matrix computation
      perf_log.pop ("Ke");

      // Now we build the element right-hand-side contribution.
      // This involves a single loop in which we integrate the
      // "forcing function" in the PDE against the test functions.
      //
      // Start logging the right-hand-side computation
      perf_log.push ("Fe");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // fxy is the forcing function for the Poisson equation.
          // In this case we set fxy to be a finite difference
          // Laplacian approximation to the (known) exact solution.
          //
          // We will use the second-order accurate FD Laplacian
          // approximation, which in 2D on a structured grid is
          //
          // u_xx + u_yy = (u(i-1,j) + u(i+1,j) +
          //                u(i,j-1) + u(i,j+1) +
          //                -4*u(i,j))/h^2
          //
          // Since the value of the forcing function depends only
          // on the location of the quadrature point (q_point[qp])
          // we will compute it here, outside of the i-loop
          const Real x = q_point[qp](0);
#if LIBMESH_DIM > 1
          const Real y = q_point[qp](1);
#else
          const Real y = 0.;
#endif
#if LIBMESH_DIM > 2
          const Real z = q_point[qp](2);
#else
          const Real z = 0.;
#endif
          const Real eps = 1.e-3;

          const Real uxx = (exact_solution(x-eps, y, z) +
                            exact_solution(x+eps, y, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uyy = (exact_solution(x, y-eps, z) +
                            exact_solution(x, y+eps, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uzz = (exact_solution(x, y, z-eps) +
                            exact_solution(x, y, z+eps) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          Real fxy;
          if (dim==1)
            {
              // In 1D, compute the rhs by differentiating the
              // exact solution twice.
              const Real pi = libMesh::pi;
              fxy = (0.25*pi*pi)*sin(.5*pi*x);
            }
          else
            {
              fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
            }

          // Add the RHS contribution
          for (std::size_t i=0; i<phi.size(); i++)
            Fe(i) += JxW[qp]*fxy*phi[i][qp];
        }

      // Stop logging the right-hand-side computation
      perf_log.pop ("Fe");

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      // Also, note that here we call heterogenously_constrain_element_matrix_and_vector
      // to impose a inhomogeneous Dirichlet boundary conditions.
      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.push ("matrix insertion");

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.pop ("matrix insertion");
    }

  // That's it.  We don't need to do anything else to the
  // PerfLog.  When it goes out of scope (at this function return)
  // it will print its log to the screen. Pretty easy, huh?
}
