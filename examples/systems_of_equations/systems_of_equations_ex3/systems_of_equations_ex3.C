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



// <h1>Systems Example 3 - Navier-Stokes with SCALAR Lagrange Multiplier</h1>
// \author David Knezevic
// \date 2010
//
// This example shows how the transient Navier-Stokes problem from
// example 13 can be solved using a scalar Lagrange multiplier
// formulation to constrain the integral of the pressure variable,
// rather than pinning the pressure at a single point.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes (EquationSystems & es,
                      const std::string & system_name);

// Function which sets Dirichlet BCs for the lid-driven cavity.
void set_lid_driven_bcs(TransientLinearImplicitSystem & system);

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // This example NaNs with the Eigen sparse linear solvers and
  // Trilinos solvers, but should work OK with either PETSc or
  // Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements in 2D.  Building these
  // higher-order elements allows us to use higher-order
  // approximation, as in example 3.
  MeshTools::Generation::build_square (mesh,
                                       20, 20,
                                       0., 1.,
                                       0., 1.,
                                       QUAD9);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Creates a transient system named "Navier-Stokes"
  TransientLinearImplicitSystem & system =
    equation_systems.add_system<TransientLinearImplicitSystem> ("Navier-Stokes");

  // Add the variables "vel_x" & "vel_y" to "Navier-Stokes".  They
  // will be approximated using second-order approximation.
  system.add_variable ("vel_x", SECOND);
  system.add_variable ("vel_y", SECOND);

  // Add the variable "p" to "Navier-Stokes". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  system.add_variable ("p", FIRST);

  // Add a scalar Lagrange multiplier to constrain the
  // pressure to have zero mean.
  system.add_variable ("alpha", FIRST, SCALAR);

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_stokes);

  // Set Dirichlet boundary conditions.
  set_lid_driven_bcs(system);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Create a performance-logging object for this example
  PerfLog perf_log("Systems Example 3");

  // Get a reference to the Stokes system to use later.
  TransientLinearImplicitSystem & navier_stokes_system =
    equation_systems.get_system<TransientLinearImplicitSystem>("Navier-Stokes");

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  const Real dt = 0.1;
  navier_stokes_system.time = 0.0;
  const unsigned int n_timesteps = 15;

  // The number of steps and the stopping criterion are also required
  // for the nonlinear iterations.
  const unsigned int n_nonlinear_steps = 15;
  const Real nonlinear_tolerance = TOLERANCE*10;

  // We also set a standard linear solver flag in the EquationSystems object
  // which controls the maxiumum number of linear solver iterations allowed.
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;

  // Tell the system of equations what the timestep is by using
  // the set_parameter function.  The matrix assembly routine can
  // then reference this parameter.
  equation_systems.parameters.set<Real> ("dt") = dt;

  // The kinematic viscosity, nu = mu/rho, units of length**2/time.
  equation_systems.parameters.set<Real> ("nu") = .007;

  // The first thing to do is to get a copy of the solution at
  // the current nonlinear iteration.  This value will be used to
  // determine if we can exit the nonlinear loop.
  UniquePtr<NumericVector<Number> >
    last_nonlinear_soln (navier_stokes_system.solution->clone());

#ifdef LIBMESH_HAVE_EXODUS_API
  // Since we are not doing adaptivity, write all solutions to a single Exodus file.
  ExodusII_IO exo_io(mesh);

  // Write out the initial condition
  exo_io.write_equation_systems ("out.e", equation_systems);
#endif

  for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
    {
      // Incremenet the time counter, set the time and the
      // time step size as parameters in the EquationSystem.
      navier_stokes_system.time += dt;

      // A pretty update message
      libMesh::out << "\n\n*** Solving time step "
                   << t_step
                   << ", time = "
                   << navier_stokes_system.time
                   << " ***"
                   << std::endl;

      // Now we need to update the solution vector from the
      // previous time step.  This is done directly through
      // the reference to the Stokes system.
      *navier_stokes_system.old_local_solution = *navier_stokes_system.current_local_solution;

      // At the beginning of each solve, reset the linear solver tolerance
      // to a "reasonable" starting value.
      const Real initial_linear_solver_tol = 1.e-6;
      equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

      // We'll set this flag when convergence is (hopefully) achieved.
      bool converged = false;

      // Now we begin the nonlinear loop
      for (unsigned int l=0; l<n_nonlinear_steps; ++l)
        {
          // Update the nonlinear solution.
          last_nonlinear_soln->zero();
          last_nonlinear_soln->add(*navier_stokes_system.solution);

          // Assemble & solve the linear system.
          perf_log.push("linear solve");
          equation_systems.get_system("Navier-Stokes").solve();
          perf_log.pop("linear solve");

          // Compute the difference between this solution and the last
          // nonlinear iterate.
          last_nonlinear_soln->add (-1., *navier_stokes_system.solution);

          // Close the vector before computing its norm
          last_nonlinear_soln->close();

          // Compute the l2 norm of the difference
          const Real norm_delta = last_nonlinear_soln->l2_norm();

          // How many iterations were required to solve the linear system?
          const unsigned int n_linear_iterations = navier_stokes_system.n_linear_iterations();

          // What was the final residual of the linear system?
          const Real final_linear_residual = navier_stokes_system.final_linear_residual();

          // If the solver did no work (sometimes -ksp_converged_reason
          // says "Linear solve converged due to CONVERGED_RTOL
          // iterations 0") but the nonlinear residual norm is above
          // the tolerance, we need to pick an even lower linear
          // solver tolerance and try again.  Note that the tolerance
          // is relative to the norm of the RHS, which for this
          // particular problem does not go to zero, since we are
          // solving for the full solution rather than the update.
          //
          // Similarly, if the solver did no work and this is the 0th
          // nonlinear step, it means that the delta between solutions
          // is being inaccurately measured as "0" since the solution
          // did not change.  Decrease the tolerance and try again.
          if (n_linear_iterations == 0 &&
              (navier_stokes_system.final_linear_residual() >= nonlinear_tolerance || l==0))
            {
              Real old_linear_solver_tolerance = equation_systems.parameters.get<Real> ("linear solver tolerance");
              equation_systems.parameters.set<Real> ("linear solver tolerance") = 1.e-3 * old_linear_solver_tolerance;
              continue;
            }

          // Print out convergence information for the linear and
          // nonlinear iterations.
          libMesh::out << "Linear solver converged at step: "
                       << n_linear_iterations
                       << ", final residual: "
                       << final_linear_residual
                       << "  Nonlinear convergence: ||u - u_old|| = "
                       << norm_delta
                       << std::endl;

          // Terminate the solution iteration if the difference between
          // this nonlinear iterate and the last is sufficiently small, AND
          // if the most recent linear system was solved to a sufficient tolerance.
          if ((norm_delta < nonlinear_tolerance) &&
              (navier_stokes_system.final_linear_residual() < nonlinear_tolerance))
            {
              libMesh::out << " Nonlinear solver converged at step "
                           << l
                           << std::endl;
              converged = true;
              break;
            }

          // Otherwise, decrease the linear system tolerance.  For the inexact Newton
          // method, the linear solver tolerance needs to decrease as we get closer to
          // the solution to ensure quadratic convergence.  The new linear solver tolerance
          // is chosen (heuristically) as the square of the previous linear system residual norm.
          //Real flr2 = final_linear_residual*final_linear_residual;
          Real new_linear_solver_tolerance = std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol);
          equation_systems.parameters.set<Real> ("linear solver tolerance") = new_linear_solver_tolerance;
        } // end nonlinear loop

      // Don't keep going if we failed to converge.
      if (!converged)
        libmesh_error_msg("Error: Newton iterations failed to converge!");

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out every nth timestep to file.
      const unsigned int write_interval = 1;

      if ((t_step+1)%write_interval == 0)
        {
          exo_io.write_timestep("out.e",
                                equation_systems,
                                t_step+1, // we're off by one since we wrote the IC and the Exodus numbering is 1-based.
                                navier_stokes_system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    } // end timestep loop.

  // All done.
  return 0;
}






// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_stokes (EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Navier-Stokes");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("Navier-Stokes");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("vel_x");
  const unsigned int v_var = navier_stokes_system.variable_number ("vel_y");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");
  const unsigned int alpha_var = navier_stokes_system.variable_number ("alpha");

  // Get the Finite Element type for "vel_x".  Note this will be
  // the same as the type for "vel_y".
  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  // Build a Finite Element object of the specified type for
  // the pressure variables.
  UniquePtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe_vel->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dphi = fe_vel->get_dphi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & psi = fe_pres->get_phi();

  // The value of the linear shape function gradients at the quadrature points
  // const std::vector<std::vector<RealGradient> > & dpsi = fe_pres->get_dphi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = navier_stokes_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke);
  DenseSubMatrix<Number> Kalpha_p(Ke), Kp_alpha(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_alpha;

  // Find out what the timestep size parameter is from the system, and
  // the value of theta for the theta method.  We use implicit Euler (theta=1)
  // for this simulation even though it is only first-order accurate in time.
  // The reason for this decision is that the second-order Crank-Nicolson
  // method is notoriously oscillatory for problems with discontinuous
  // initial data such as the lid-driven cavity.  Therefore,
  // we sacrifice accuracy in time for stability, but since the solution
  // reaches steady state relatively quickly we can afford to take small
  // timesteps.  If you monitor the initial nonlinear residual for this
  // simulation, you should see that it is monotonically decreasing in time.
  const Real dt    = es.parameters.get<Real>("dt");
  const Real theta = 1.;

  // The kinematic viscosity, multiplies the "viscous" terms.
  const Real nu = es.parameters.get<Real>("nu");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
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
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_alpha, alpha_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Reposition the submatrices...  The idea is this:
      //
      //         -           -          -  -
      //        | Kuu Kuv Kup |        | Fu |
      //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
      //        | Kpu Kpv Kpp |        | Fp |
      //         -           -          -  -
      //
      // The DenseSubMatrix.repostition () member takes the
      // (row_offset, column_offset, row_size, column_size).
      //
      // Similarly, the DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      // Also, add a row and a column to constrain the pressure
      Kp_alpha.reposition (p_var*n_u_dofs, p_var*n_u_dofs+n_p_dofs, n_p_dofs, 1);
      Kalpha_p.reposition (p_var*n_u_dofs+n_p_dofs, p_var*n_u_dofs, 1, n_p_dofs);


      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This must be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.
          Number u = 0., u_old = 0.;
          Number v = 0., v_old = 0.;
          Number p_old = 0.;
          Gradient grad_u, grad_u_old;
          Gradient grad_v, grad_v_old;

          // Compute the velocity & its gradient from the previous timestep
          // and the old Newton iterate.
          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              // From the old timestep:
              u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              grad_u_old.add_scaled (dphi[l][qp], navier_stokes_system.old_solution (dof_indices_u[l]));
              grad_v_old.add_scaled (dphi[l][qp], navier_stokes_system.old_solution (dof_indices_v[l]));

              // From the previous Newton iterate:
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              grad_u.add_scaled (dphi[l][qp], navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp], navier_stokes_system.current_solution (dof_indices_v[l]));
            }

          // Compute the old pressure value at this quadrature point.
          for (unsigned int l=0; l<n_p_dofs; l++)
            p_old += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);

          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (u_old, v_old);
          const NumberVectorValue U     (u,     v);
          const Number u_x = grad_u(0);
          const Number u_y = grad_u(1);
          const Number v_x = grad_v(0);
          const Number v_y = grad_v(1);

          // First, an i-loop over the velocity degrees of freedom.
          // We know that n_u_dofs == n_v_dofs so we can compute contributions
          // for both at the same time.
          for (unsigned int i=0; i<n_u_dofs; i++)
            {
              Fu(i) += JxW[qp]*(u_old*phi[i][qp] -                            // mass-matrix term
                                (1.-theta)*dt*(U_old*grad_u_old)*phi[i][qp] + // convection term
                                (1.-theta)*dt*p_old*dphi[i][qp](0)  -         // pressure term on rhs
                                (1.-theta)*dt*nu*(grad_u_old*dphi[i][qp]) +   // diffusion term on rhs
                                theta*dt*(U*grad_u)*phi[i][qp]);              // Newton term


              Fv(i) += JxW[qp]*(v_old*phi[i][qp] -                             // mass-matrix term
                                (1.-theta)*dt*(U_old*grad_v_old)*phi[i][qp] +  // convection term
                                (1.-theta)*dt*p_old*dphi[i][qp](1) -           // pressure term on rhs
                                (1.-theta)*dt*nu*(grad_v_old*dphi[i][qp]) +    // diffusion term on rhs
                                theta*dt*(U*grad_v)*phi[i][qp]);               // Newton term


              // Note that the Fp block is identically zero unless we are using
              // some kind of artificial compressibility scheme...

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j<n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                 // mass matrix term
                                       theta*dt*nu*(dphi[i][qp]*dphi[j][qp]) + // diffusion term
                                       theta*dt*(U*dphi[j][qp])*phi[i][qp] +   // convection term
                                       theta*dt*u_x*phi[i][qp]*phi[j][qp]);    // Newton term

                  Kuv(i,j) += JxW[qp]*theta*dt*u_y*phi[i][qp]*phi[j][qp];      // Newton term

                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                 // mass matrix term
                                       theta*dt*nu*(dphi[i][qp]*dphi[j][qp]) + // diffusion term
                                       theta*dt*(U*dphi[j][qp])*phi[i][qp] +   // convection term
                                       theta*dt*v_y*phi[i][qp]*phi[j][qp]);    // Newton term

                  Kvu(i,j) += JxW[qp]*theta*dt*v_x*phi[i][qp]*phi[j][qp];      // Newton term
                }

              // Matrix contributions for the up and vp couplings.
              for (unsigned int j=0; j<n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](0));
                  Kvp(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](1));
                }
            }

          // Now an i-loop over the pressure degrees of freedom.  This code computes
          // the matrix entries due to the continuity equation.  Note: To maintain a
          // symmetric matrix, we may (or may not) multiply the continuity equation by
          // negative one.  Here we do not.
          for (unsigned int i=0; i<n_p_dofs; i++)
            {
              Kp_alpha(i,0) += JxW[qp]*psi[i][qp];
              Kalpha_p(0,i) += JxW[qp]*psi[i][qp];
              for (unsigned int j=0; j<n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                  Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                }
            }
        } // end of the quadrature point qp-loop

      // Since we're using heterogeneous DirichletBoundary objects for
      // the boundary conditions, we need to call a specific function
      // to constrain the element stiffness matrix.
      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  // We can set the mean of the pressure by setting Falpha.  Typically
  // a value of zero is chosen, but the value should be arbitrary.
  navier_stokes_system.rhs->add(navier_stokes_system.rhs->size()-1, 10.);
}



void set_lid_driven_bcs(TransientLinearImplicitSystem & system)
{
  unsigned short int
    u_var = system.variable_number("vel_x"),
    v_var = system.variable_number("vel_y");

  // Get a convenient reference to the System's DofMap
  DofMap & dof_map = system.get_dof_map();

  {
    // u=v=0 on bottom, left, right
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(0);
    boundary_ids.insert(1);
    boundary_ids.insert(3);

    std::vector<unsigned int> variables;
    variables.push_back(u_var);
    variables.push_back(v_var);

    dof_map.add_dirichlet_boundary(DirichletBoundary(boundary_ids,
                                                     variables,
                                                     ZeroFunction<Number>()));
  }
  {
    // u=1 on top
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(2);

    std::vector<unsigned int> variables;
    variables.push_back(u_var);

    dof_map.add_dirichlet_boundary(DirichletBoundary(boundary_ids,
                                                     variables,
                                                     ConstFunction<Number>(1.)));
  }
  {
    // v=0 on top
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(2);

    std::vector<unsigned int> variables;
    variables.push_back(v_var);

    dof_map.add_dirichlet_boundary(DirichletBoundary(boundary_ids,
                                                     variables,
                                                     ZeroFunction<Number>()));
  }
}
