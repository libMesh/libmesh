/* $Id: ex18.C,v 1.5 2006-06-29 19:38:43 roystgnr Exp $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

 // <h1>Example 18 - Unsteady Navier-Stokes Equations with DiffSystem</h1>
 //
 // This example shows how the transient nonlinear problem from
 // example 13 can be solved using the DiffSystem class framework to
 // simplify the user-implemented equations.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files
#include "equation_systems.h"
#include "error_vector.h"
#include "fe_base.h"
#include "gmv_io.h"
#include "kelly_error_estimator.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_refinement.h"
#include "quadrature.h"
#include "uniform_refinement_estimator.h"

// Some (older) compilers do not offer full stream 
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// DiffSystem framework files
#include "diff_solver.h"
#include "fem_system.h"
#include "euler_solver.h"
#include "steady_solver.h"

// The global FEM error tolerance at each timestep
// Make this nonzero to solve to a specified tolerance
// This will probably break with KellyErrorIndicator
// const Real global_tolerance = 1.e-3;
const Real global_tolerance = 0.;

// The desired number of active mesh elements
// Make this nonzero to solve to a specified mesh size
const unsigned int nelem_target = 1000;

// Solve a transient instead of a steady problem?
const bool transient = true;

// The interval between our timesteps
const Real deltat = 0.005;

// And the number of timesteps to take
unsigned int n_timesteps = 15;

// Write out every nth timestep to file.
const unsigned int write_interval = 1;

// The coarse grid size from which to start adaptivity
const unsigned int coarsegridsize = 1;

// The maximum number of adaptive steps per timestep
const unsigned int n_adaptivesteps = 10;

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class NavierSystem : public FEMSystem
{
public:
  // Constructor
  NavierSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number)
  : FEMSystem(es, name, number), Reynolds(0.) {}

  // System initialization
  virtual void init_data ();

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian);

  // Constraint parts
  virtual bool element_constraint (bool request_jacobian);
  virtual bool side_constraint (bool request_jacobian);

  // Finite elements for the velocity and pressure on element interiors
  FEBase *fe_velocity, *fe_pressure;

  // Finite element for the velocity on element sides
  FEBase *fe_side_vel;

  // The Reynolds number to solve for
  Real Reynolds;
};


// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  libMesh::init (argc, argv);
  {    
    // Create a two-dimensional mesh.
    Mesh mesh (2);
    
    // Use the MeshTools::Generation mesh generator to create a uniform
    // grid on the square [-1,1]^D.  We instruct the mesh generator
    // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
    // elements in 3D.  Building these higher-order elements allows
    // us to use higher-order approximation, as in example 3.
    MeshTools::Generation::build_square (mesh,
                                         coarsegridsize,
                                         coarsegridsize,
                                         0., 1.,
                                         0., 1.,
                                         QUAD9);
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
    // Declare the system "Navier-Stokes" and its variables.
    NavierSystem & stokes_system = 
      equation_systems.add_system<NavierSystem> ("Navier-Stokes");

    // Solve this as a time-dependent or steady system
    if (transient)
      stokes_system.time_solver =
        AutoPtr<TimeSolver>(new EulerSolver(stokes_system));
    else
      {
        stokes_system.time_solver =
          AutoPtr<TimeSolver>(new SteadySolver(stokes_system));
        assert(n_timesteps == 1);
      }

    // Initialize the system
    equation_systems.init ();

    // Set the time stepping options
    stokes_system.deltat = deltat;

    // And the nonlinear solver options
    DiffSolver &solver = *stokes_system.time_solver->diff_solver;
    // solver.quiet = false;
    solver.max_nonlinear_iterations = 15;
    solver.relative_step_tolerance = 1.e-7;
    solver.relative_residual_tolerance = 1.e-8;

    // And the linear solver options
    solver.max_linear_iterations = 1000;
    solver.initial_linear_tolerance = 1.e-3;

    // Print information about the system to the screen.
    equation_systems.print_info();

    MeshRefinement mesh_refinement(mesh);

    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
      {
        // A pretty update message
        std::cout << " Solving time step " << t_step << ", time = "
                  << stokes_system.time << std::endl;

        // Adaptively solve the timestep
        unsigned int a_step = 0;
        for (; a_step != n_adaptivesteps; ++a_step)
          {
            stokes_system.solve();

            ErrorVector error;

            AutoPtr<ErrorEstimator> error_estimator;

            // To solve to a tolerance in this problem we
            // need a better estimator than Kelly
            if (global_tolerance != 0.)
              {
                // We can't adapt to both a tolerance and a mesh
                // size at once
                assert (nelem_target == 0);

                UniformRefinementEstimator *u =
                  new UniformRefinementEstimator;

                // The lid-driven cavity problem isn't in H1, so
                // lets estimate H0 (i.e. L2) error
                u->sobolev_order() = 0;

                error_estimator.reset(u);
              }
            else
              {
		// If we aren't adapting to a tolerance we need a
                // target mesh size
                assert (nelem_target > 0);

                error_estimator.reset(new KellyErrorEstimator);
              }

            // Calculate error based on u and v but not p
            error_estimator->component_scale.push_back(1.0); // u
            error_estimator->component_scale.push_back(1.0); // v
            error_estimator->component_scale.push_back(0.0); // p

            error_estimator->estimate_error(stokes_system, error);

            // Print out status at each adaptive step.
            Real global_error = error.l2_norm();
            std::cerr << "adaptive step " << a_step << ": ";
            if (global_tolerance != 0.)
              std::cerr << "global_error = " << global_error
                        << " with ";
            std::cerr << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
                      << " active dofs." << std::endl;
            if (global_tolerance != 0.)
              std::cerr << "worst element error = " << error.maximum()
                        << ", mean = " << error.mean() << std::endl;

            if (global_tolerance != 0.)
              {
                // If we've reached our desired tolerance, we
                // don't need any more adaptive steps
                if (global_error < global_tolerance)
                  break;
                mesh_refinement.flag_elements_by_error_tolerance
                  (error, global_tolerance);
              }
            else
              {
		// If flag_elements_to_nelem_target returns true, this
                // should be our last adaptive step.
                if (mesh_refinement.flag_elements_to_nelem_target
                      (error, nelem_target))
                  {
                    mesh_refinement.refine_and_coarsen_elements();
                    equation_systems.reinit();
                    break;
                  }
              }

            // Carry out the adaptive mesh refinement/coarsening
            mesh_refinement.refine_and_coarsen_elements();
            equation_systems.reinit();
          }
        // Do one last solve if necessary
        if (a_step == n_adaptivesteps)
          {
            stokes_system.solve();
          }

        // Advance to the next timestep in a transient problem
        stokes_system.time_solver->advance_timestep();

        // Write out this timestep if we're requested to
        if ((t_step+1)%write_interval == 0)
          {
            OStringStream file_name;

            // We write the file name in the gmv auto-read format.
            file_name << "out.gmv.";
            OSSRealzeroright(file_name,3,0, t_step + 1);

            GMVIO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
          }
      }
  }
  
  // All done.  
  return libMesh::close ();
}



void NavierSystem::init_data ()
{
  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  this->add_variable ("u", SECOND);
  this->add_variable ("v", SECOND);

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  this->add_variable ("p", FIRST);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Tell the system to march u and v forward in time, but 
  // leave p as a constraint only
  this->time_evolving(0);
  this->time_evolving(1);

  // Get references to the finite elements we need
  fe_velocity = element_fe[this->variable_type(0)];
  fe_pressure = element_fe[this->variable_type(2)];
  fe_side_vel = side_fe[this->variable_type(0)];

  // To enable FE optimizations, we should prerequest all the data
  // we will need to build the linear system.
  fe_velocity->get_JxW();
  fe_velocity->get_phi();
  fe_velocity->get_dphi();

  fe_pressure->get_phi();

  fe_side_vel->get_JxW();
  fe_side_vel->get_phi();
  fe_side_vel->get_xyz();

  // Useful debugging options
  // this->verify_analytic_jacobians = 1e-6;
  // this->print_jacobians = true;
}


bool NavierSystem::element_time_derivative (bool request_jacobian)
{
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = fe_velocity->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> > &phi = fe_velocity->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> > &dphi =
    fe_velocity->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> > &psi = fe_pressure->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = dof_indices_var[0].size(); 
  assert (n_u_dofs == dof_indices_var[1].size()); 
  const unsigned int n_p_dofs = dof_indices_var[2].size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &Kuu = *elem_subjacobians[0][0];
  DenseSubMatrix<Number> &Kvv = *elem_subjacobians[1][1];
  DenseSubMatrix<Number> &Kuv = *elem_subjacobians[0][1];
  DenseSubMatrix<Number> &Kvu = *elem_subjacobians[1][0];
  DenseSubMatrix<Number> &Kup = *elem_subjacobians[0][2];
  DenseSubMatrix<Number> &Kvp = *elem_subjacobians[1][2];
  DenseSubVector<Number> &Fu = *elem_subresiduals[0];
  DenseSubVector<Number> &Fv = *elem_subresiduals[1];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      Number u = interior_value(0, qp),
             v = interior_value(1, qp),
             p = interior_value(2, qp);
      Gradient grad_u = interior_gradient(0, qp),
               grad_v = interior_gradient(1, qp);

      // Definitions for convenience.  It is sometimes simpler to do a
      // dot product if you have the full vector at your disposal.
      const NumberVectorValue U     (u,     v);
      const Number  u_x = grad_u(0);
      const Number  u_y = grad_u(1);
      const Number  v_x = grad_v(0);
      const Number  v_y = grad_v(1);
          
      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] *
                   (-Reynolds*(U*grad_u)*phi[i][qp] + // convection term
                    p*dphi[i][qp](0) -                // pressure term
                    (grad_u*dphi[i][qp]));            // diffusion term
            
          Fv(i) += JxW[qp] *
                   (-Reynolds*(U*grad_v)*phi[i][qp] + // convection term
                    p*dphi[i][qp](1) -                // pressure term
                    (grad_v*dphi[i][qp]));            // diffusion term

          // Note that the Fp block is identically zero unless we are using
          // some kind of artificial compressibility scheme...

          // Matrix contributions for the uu and vv couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW[qp] *
 /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
 /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
 /* Newton term     */       Reynolds*u_x*phi[i][qp]*phi[j][qp]);

                Kuv(i,j) += JxW[qp] *
 /* Newton term     */      -Reynolds*u_y*phi[i][qp]*phi[j][qp];

                Kvv(i,j) += JxW[qp] *
 /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
 /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
 /* Newton term     */       Reynolds*v_y*phi[i][qp]*phi[j][qp]);

                Kvu(i,j) += JxW[qp] * 
 /* Newton term     */      -Reynolds*v_x*phi[i][qp]*phi[j][qp];
            }

          // Matrix contributions for the up and vp couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_p_dofs; j++)
              {
                Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
              }
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}



bool NavierSystem::element_constraint (bool request_jacobian)
{
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> &JxW = fe_velocity->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> > &dphi =
    fe_velocity->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> > &psi = fe_pressure->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = dof_indices_var[0].size();
  const unsigned int n_p_dofs = dof_indices_var[2].size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &Kpu = *elem_subjacobians[2][0];
  DenseSubMatrix<Number> &Kpv = *elem_subjacobians[2][1];
  DenseSubVector<Number> &Fp = *elem_subresiduals[2];

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      Gradient grad_u = interior_gradient(0, qp),
               grad_v = interior_gradient(1, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * psi[i][qp] *
                   (grad_u(0) + grad_v(1));

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
              }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

      

bool NavierSystem::side_constraint (bool request_jacobian)
{
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = fe_side_vel->get_JxW();

  // The velocity shape functions at side quadrature points.
  const std::vector<std::vector<Real> > &phi_side =
    fe_side_vel->get_phi();

  // The XYZ locations (in physical space) of the
  // quadrature points on the side.  This is where
  // we will interpolate the boundary value function.
  const std::vector<Point > &qside_point = fe_side_vel->get_xyz();

  // The number of local degrees of freedom in u and v
  const unsigned int n_u_dofs = dof_indices_var[0].size(); 

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &Kuu = *elem_subjacobians[0][0];
  DenseSubMatrix<Number> &Kvv = *elem_subjacobians[1][1];

  DenseSubVector<Number> &Fu = *elem_subresiduals[0];
  DenseSubVector<Number> &Fv = *elem_subresiduals[1];

  // For this example we will use Dirichlet velocity boundary
  // conditions imposed at each timestep via the penalty method.

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const Real penalty = 1.e10;

  unsigned int n_sidepoints = side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      Number u = side_value(0, qp),
             v = side_value(1, qp);

      // The location of the current boundary quadrature point
      // const Real xf = qside_point[qp](0);
      const Real yf = qside_point[qp](1);

      // Set u = 1 on the top boundary, 0 everywhere else
      const Real u_value = (yf > 1.0 - TOLERANCE) ? 1. : 0.;
              
      // Set v = 0 everywhere
      const Real v_value = 0.;

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * penalty *
                   (u - u_value) * phi_side[i][qp];
          Fv(i) += JxW_side[qp] * penalty *
                   (v - v_value) * phi_side[i][qp];

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW_side[qp] * penalty *
                            phi_side[i][qp] * phi_side[j][qp];
                Kvv(i,j) += JxW_side[qp] * penalty *
                            phi_side[i][qp] * phi_side[j][qp];
              }
        }
    }

  return request_jacobian;
}
