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



// <h1>Miscellaneous Example 3 - 2D Laplace-Young Problem Using Nonlinear Solvers</h1>
// \author Derek Gaston
// \date 2008
//
// This example shows how to use the NonlinearImplicitSystem class to
// solve nonlinear problems in libMesh.  The NonlinearImplicitSystem
// class employs the user's ComputeResidual and ComputeJacobian
// objects during the solve.  In this particular example, the
// LaplaceYoung object provides this functionality.
//
// You can turn on preconditioning of the matrix-free system using the
// jacobian by passing "-pre" on the command line.  Currently, this
// feature only works with Petsc, so this isn't used by "make run".
//
// This example also runs with the experimental Trilinos NOX solvers
// by specifying the --use-trilinos command line argument.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

// Necessary for programmatically setting petsc options
#ifdef LIBMESH_HAVE_PETSC
#include <petsc.h>
#include "libmesh/petsc_macro.h"
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;

/**
 * A class which provides the residual and jacobian assembly
 * functions for the Laplace-Young system of equations.
 */
class LaplaceYoung :
  public NonlinearImplicitSystem::ComputeResidual,
  public NonlinearImplicitSystem::ComputeJacobian,
  public NonlinearImplicitSystem::ComputePostCheck
{
public:
  LaplaceYoung() :
    _kappa(1.),
    _sigma(0.2),
    _gamma(1.0)
  {}

  /**
   * Function which computes the residual.
   */
  virtual void residual (const NumericVector<Number> & X,
                         NumericVector<Number> & R,
                         NonlinearImplicitSystem & S);

  /**
   * Function which computes the jacobian.
   */
  virtual void jacobian (const NumericVector<Number> & X,
                         SparseMatrix<Number> & J,
                         NonlinearImplicitSystem & S);

  /**
   * Function which performs a postcheck on the solution.  In this
   * example, we take a "damped" Newton step defined by:
   *
   * u_new = u_old + gamma*delta_u
   *
   * where delta_u is the search direction, and 0 < gamma <= 1 is a
   * damping parameter.  gamma=1 corresponds to a full Newton step.
   *
   * This is really for demonstration purposes only, as it just
   * degrades the rate of nonlinear convergence in this particular
   * example.
   */
  virtual void postcheck (const NumericVector<Number> & old_soln,
                          NumericVector<Number> & search_direction,
                          NumericVector<Number> & new_soln,
                          bool & changed_search_direction,
                          bool & changed_new_soln,
                          NonlinearImplicitSystem & S);

private:
  Real _kappa;
  Real _sigma;

  // Damping factor used for the solve postcheck
  Real _gamma;
};



// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires a NonlinearSolver.
#if !defined(LIBMESH_HAVE_PETSC) && (!defined(LIBMESH_TRILINOS_HAVE_NOX) || !defined(LIBMESH_TRILINOS_HAVE_EPETRA))
  libmesh_example_requires(false, "--enable-petsc or --enable-trilinos");
#endif

  if (libMesh::on_command_line ("--use-eigen"))
    {
      libMesh::err << "This example requires a NonlinearSolver, and therefore does not "
                   << "support --use-eigen on the command line."
                   << std::endl;
      return 0;
    }

#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Create a GetPot object to parse the command line
  GetPot command_line (argc, argv);

  // Check for proper calling arguments.
  if (argc < 3)
    {
      // This handy function will print the file name, line number,
      // specified message, and then throw an exception.
      libmesh_error_msg("Usage:\n" << "\t " << argv[0] << " -r 2");
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


  // Read number of refinements
  int nr = 2;
  if (command_line.search(1, "-r"))
    nr = command_line.next(nr);

  // Read FE order from command line
  std::string order = "FIRST";
  if (command_line.search(2, "-Order", "-o"))
    order = command_line.next(order);

  // Read FE Family from command line
  std::string family = "LAGRANGE";
  if (command_line.search(2, "-FEFamily", "-f"))
    family = command_line.next(family);

  // Cannot use dicontinuous basis.
  if ((family == "MONOMIAL") || (family == "XYZ"))
    libmesh_error_msg("This example requires a C^0 (or higher) FE basis.");

  if (command_line.search(1, "-pre"))
    {
#ifdef LIBMESH_HAVE_PETSC
      //Use the jacobian for preconditioning.
#  if PETSC_VERSION_LESS_THAN(3,7,0)
      PetscOptionsSetValue("-snes_mf_operator", PETSC_NULL);
#  else
      PetscOptionsSetValue(PETSC_NULL, "-snes_mf_operator", PETSC_NULL);
#  endif
#else
      libMesh::err << "Must be using PETSc to use jacobian based preconditioning" << std::endl;

      //returning zero so that "make run" won't fail if we ever enable this capability there.
      return 0;
#endif //LIBMESH_HAVE_PETSC
    }

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden by the file,
  // distributed across the default MPI communicator.
  Mesh mesh(init.comm());

  mesh.read ("lshaped.xda");

  if (order != "FIRST")
    mesh.all_second_order();

  MeshRefinement(mesh).uniformly_refine(nr);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.

  // Creates a system named "Laplace-Young"
  NonlinearImplicitSystem & system =
    equation_systems.add_system<NonlinearImplicitSystem> ("Laplace-Young");

  // Here we specify the tolerance for the nonlinear solver and
  // the maximum of nonlinear iterations.
  equation_systems.parameters.set<Real>         ("nonlinear solver tolerance")          = 1.e-12;
  equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations") = 50;


  // Adds the variable "u" to "Laplace-Young".  "u"
  // will be approximated using second-order approximation.
  system.add_variable("u",
                      Utility::string_to_enum<Order>   (order),
                      Utility::string_to_enum<FEFamily>(family));

  // Consruct object which provides the residual and jacobian
  // computations and tell the solver to use it.
  LaplaceYoung laplace_young;
  system.nonlinear_solver->residual_object = &laplace_young;
  system.nonlinear_solver->jacobian_object = &laplace_young;
  system.nonlinear_solver->postcheck_object = &laplace_young;

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Solve the system "Laplace-Young", print the number of iterations
  // and final residual
  equation_systems.get_system("Laplace-Young").solve();

  // Print out final convergence information.  This duplicates some
  // output from during the solve itself, but demonstrates another way
  // to get this information after the solve is complete.
  libMesh::out << "Laplace-Young system solved at nonlinear iteration "
               << system.n_nonlinear_iterations()
               << " , final nonlinear residual norm: "
               << system.final_nonlinear_residual()
               << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API
  // After solving the system write the solution
  ExodusII_IO (mesh).write_equation_systems ("out.e",
                                             equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}



// Residual assembly function for the Laplace-Young system
void LaplaceYoung::residual (const NumericVector<Number> & soln,
                             NumericVector<Number> & residual,
                             NonlinearImplicitSystem & sys)
{
  EquationSystems & es = sys.get_equation_systems();

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  libmesh_assert_equal_to (dim, 2);

  // Get a reference to the NonlinearImplicitSystem we are solving
  NonlinearImplicitSystem & system =
    es.get_system<NonlinearImplicitSystem>("Laplace-Young");

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

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Define data structures to contain the resdual contributions
  DenseVector<Number> Re;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the active elements in the mesh which
  // are local to this processor.
  // We will compute the element residual.
  residual.zero();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
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

      // We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Re.resize (dof_indices.size());

      // Now we will build the residual. This involves
      // the construction of the matrix K and multiplication of it
      // with the current solution x. We rearrange this into two loops:
      // In the first, we calculate only the contribution of
      // K_ij*x_j which is independent of the row i. In the second loops,
      // we multiply with the row-dependent part and add it to the element
      // residual.

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number u = 0;
          Gradient grad_u;

          for (std::size_t j=0; j<phi.size(); j++)
            {
              u      += phi[j][qp]*soln(dof_indices[j]);
              grad_u += dphi[j][qp]*soln(dof_indices[j]);
            }

          const Number K = 1./std::sqrt(1. + grad_u*grad_u);

          for (std::size_t i=0; i<phi.size(); i++)
            Re(i) += JxW[qp]*(
                              K*(dphi[i][qp]*grad_u) +
                              _kappa*phi[i][qp]*u
                              );
        }

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.

      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int side=0; side<elem->n_sides(); side++)
        if (elem->neighbor_ptr(side) == libmesh_nullptr)
          {
            // The value of the shape functions at the quadrature
            // points.
            const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();

            // The Jacobian * Quadrature Weight at the quadrature
            // points on the face.
            const std::vector<Real> & JxW_face = fe_face->get_JxW();

            // Compute the shape function values on the element face.
            fe_face->reinit(elem, side);

            // Loop over the face quadrature points for integration.
            for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                // This is the right-hand-side contribution (f),
                // which has to be subtracted from the current residual
                for (std::size_t i=0; i<phi_face.size(); i++)
                  Re(i) -= JxW_face[qp]*_sigma*phi_face[i][qp];
              }
          }

      dof_map.constrain_element_vector (Re, dof_indices);
      residual.add_vector (Re, dof_indices);
    }

  // That's it.
}



// Jacobian assembly function for the Laplace-Young system
void LaplaceYoung::jacobian (const NumericVector<Number> & soln,
                             SparseMatrix<Number> & jacobian,
                             NonlinearImplicitSystem & sys)
{
  // Get a reference to the equation system.
  EquationSystems & es = sys.get_equation_systems();

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the NonlinearImplicitSystem we are solving
  NonlinearImplicitSystem & system =
    es.get_system<NonlinearImplicitSystem>("Laplace-Young");

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

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Define data structures to contain the Jacobian element matrix.
  // Following basic finite element terminology we will denote these
  // "Ke". More detail is in example 3.
  DenseMatrix<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the active elements in the mesh which
  // are local to this processor.
  // We will compute the element Jacobian contribution.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
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

      // Zero the element Jacobian before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      // Now we will build the element Jacobian.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j). Note that the Jacobian depends
      // on the current solution x, which we access using the soln
      // vector.
      //
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Gradient grad_u;

          for (std::size_t i=0; i<phi.size(); i++)
            grad_u += dphi[i][qp]*soln(dof_indices[i]);

          const Number
            sa = 1. + grad_u*grad_u,
            K  = 1. / std::sqrt(sa),
            dK = -K / sa;

          for (std::size_t i=0; i<phi.size(); i++)
            for (std::size_t j=0; j<phi.size(); j++)
              Ke(i,j) += JxW[qp]*(
                                  K * (dphi[i][qp]*dphi[j][qp]) +
                                  dK * (grad_u*dphi[j][qp]) * (grad_u*dphi[i][qp]) +
                                  _kappa * phi[i][qp] * phi[j][qp]
                                  );
        }

      dof_map.constrain_element_matrix (Ke, dof_indices);

      // Add the element matrix to the system Jacobian.
      jacobian.add_matrix (Ke, dof_indices);
    }

  // That's it.
}



// Jacobian assembly function for the Laplace-Young system
void LaplaceYoung::postcheck (const NumericVector<Number> & old_soln,
                              NumericVector<Number> & search_direction,
                              NumericVector<Number> & new_soln,
                              bool & /*changed_search_direction*/,
                              bool & changed_new_soln,
                              NonlinearImplicitSystem & /*S*/)
{
  // Back up along the search direction by some amount.  Since Newton
  // already works well for this problem, the only affect of this is
  // to degrade the rate of convergence.
  //
  // The minus sign is due to the sign of the "search_direction"
  // vector which comes in from the Newton solve. The RHS of the
  // nonlinear system, i.e. the residual, is generally not multiplied
  // by -1, so the solution vector, i.e. the search_direction, has a
  // factor -1.
  if (_gamma != 1.0)
    {
      new_soln = old_soln;
      new_soln.add(-_gamma, search_direction);
      changed_new_soln = true;
    }
  else
    changed_new_soln = false;
}
