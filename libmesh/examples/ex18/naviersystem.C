/* $Id: naviersystem.C,v 1.1 2006-06-29 20:49:42 roystgnr Exp $ */

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

#include "naviersystem.h"

#include "fe_base.h"
#include "quadrature.h"


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
