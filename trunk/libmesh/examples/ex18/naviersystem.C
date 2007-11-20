/* $Id$ */

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

#include "getpot.h"

#include "naviersystem.h"

#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"


void NavierSystem::init_data ()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  this->add_variable ("u", SECOND);
  u_var = 0;
  this->add_variable ("v", SECOND);
  v_var = 1;

  if (dim == 3)
    {
      this->add_variable ("w", SECOND);
      w_var = 2;
      p_var = 3;
    }
  else
    {
      w_var = u_var;
      p_var = 2;
    }

  this->add_variable ("p", FIRST);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Tell the system to march velocity forward in time, but 
  // leave p as a constraint only
  this->time_evolving(u_var);
  this->time_evolving(v_var);
  if (dim == 3)
    this->time_evolving(w_var);

  // Get references to the finite elements we need
  fe_velocity = element_fe[this->variable_type(u_var)];
  fe_pressure = element_fe[this->variable_type(p_var)];
  fe_side_vel = side_fe[this->variable_type(u_var)];

  // To enable FE optimizations, we should prerequest all the data
  // we will need to build the linear system.
  fe_velocity->get_JxW();
  fe_velocity->get_phi();
  fe_velocity->get_dphi();
  fe_velocity->get_xyz();
  
  fe_pressure->get_phi();

  fe_side_vel->get_JxW();
  fe_side_vel->get_phi();

  // Check the input file for Reynolds number
  GetPot infile("navier.in");
  Reynolds = infile("Reynolds", 1.);
  application = infile("application", 0);

  // Useful debugging options
  // Set verify_analytic_jacobians to 1e-6 to use
  this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0);
  this->print_jacobians = infile("print_jacobians", false);
  this->print_element_jacobians = infile("print_element_jacobians", false);
}


bool NavierSystem::element_time_derivative (bool request_jacobian)
{
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = fe_velocity->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = fe_velocity->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    fe_velocity->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pressure->get_phi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = fe_velocity->get_xyz();
 
  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = dof_indices_var[p_var].size();
  const unsigned int n_u_dofs = dof_indices_var[u_var].size(); 
  assert (n_u_dofs == dof_indices_var[v_var].size()); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *elem_subjacobians[w_var][w_var];
  DenseSubMatrix<Number> &Kuv = *elem_subjacobians[u_var][v_var];
  DenseSubMatrix<Number> &Kuw = *elem_subjacobians[u_var][w_var];
  DenseSubMatrix<Number> &Kvu = *elem_subjacobians[v_var][u_var];
  DenseSubMatrix<Number> &Kvw = *elem_subjacobians[v_var][w_var];
  DenseSubMatrix<Number> &Kwu = *elem_subjacobians[w_var][u_var];
  DenseSubMatrix<Number> &Kwv = *elem_subjacobians[w_var][v_var];
  DenseSubMatrix<Number> &Kup = *elem_subjacobians[u_var][p_var];
  DenseSubMatrix<Number> &Kvp = *elem_subjacobians[v_var][p_var];
  DenseSubMatrix<Number> &Kwp = *elem_subjacobians[w_var][p_var];
  DenseSubVector<Number> &Fu = *elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *elem_subresiduals[w_var];
      
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
      Number p = interior_value(p_var, qp),
             u = interior_value(u_var, qp),
             v = interior_value(v_var, qp),
             w = interior_value(w_var, qp);
      Gradient grad_u = interior_gradient(u_var, qp),
               grad_v = interior_gradient(v_var, qp),
               grad_w = interior_gradient(w_var, qp);

      // Definitions for convenience.  It is sometimes simpler to do a
      // dot product if you have the full vector at your disposal.
      NumberVectorValue U     (u,     v);
      if (dim == 3)
        U(2) = w;
      const Number  u_x = grad_u(0);
      const Number  u_y = grad_u(1);
      const Number  u_z = (dim == 3)?grad_u(2):0;
      const Number  v_x = grad_v(0);
      const Number  v_y = grad_v(1);
      const Number  v_z = (dim == 3)?grad_v(2):0;
      const Number  w_x = (dim == 3)?grad_w(0):0;
      const Number  w_y = (dim == 3)?grad_w(1):0;
      const Number  w_z = (dim == 3)?grad_w(2):0;

      // Value of the forcing function at this quadrature point
      Point f = this->forcing(qpoint[qp]);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] *
                   (-Reynolds*(U*grad_u)*phi[i][qp] + // convection term
                    p*dphi[i][qp](0) -                // pressure term
		    (grad_u*dphi[i][qp]) +            // diffusion term
		    f(0)*phi[i][qp]                   // forcing function
		    );

            
          Fv(i) += JxW[qp] *
                   (-Reynolds*(U*grad_v)*phi[i][qp] + // convection term
                    p*dphi[i][qp](1) -                // pressure term
		    (grad_v*dphi[i][qp]) +            // diffusion term
		    f(1)*phi[i][qp]                   // forcing function
		    );


          if (dim == 3)
          Fw(i) += JxW[qp] *
                   (-Reynolds*(U*grad_w)*phi[i][qp] + // convection term
                    p*dphi[i][qp](2) -                // pressure term
		    (grad_w*dphi[i][qp]) +            // diffusion term
		    f(2)*phi[i][qp]                   // forcing function
		    );


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
                if (dim == 3)
                  {
                    Kww(i,j) += JxW[qp] *
 /* convection term */          (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
 /* diffusion term  */           (dphi[i][qp]*dphi[j][qp]) -
 /* Newton term     */           Reynolds*w_z*phi[i][qp]*phi[j][qp]);
                    Kuw(i,j) += JxW[qp] *
 /* Newton term     */      -Reynolds*u_z*phi[i][qp]*phi[j][qp];
                    Kvw(i,j) += JxW[qp] *
 /* Newton term     */      -Reynolds*v_z*phi[i][qp]*phi[j][qp];
                    Kwu(i,j) += JxW[qp] *
 /* Newton term     */      -Reynolds*w_x*phi[i][qp]*phi[j][qp];
                    Kwv(i,j) += JxW[qp] *
 /* Newton term     */      -Reynolds*w_y*phi[i][qp]*phi[j][qp];
                  }
              }

          // Matrix contributions for the up and vp couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_p_dofs; j++)
              {
                Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
                if (dim == 3)
                  Kwp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);
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
  const std::vector<std::vector<RealGradient> >& dphi =
    fe_velocity->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pressure->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = dof_indices_var[u_var].size();
  const unsigned int n_p_dofs = dof_indices_var[p_var].size();

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kpu = *elem_subjacobians[p_var][u_var];
  DenseSubMatrix<Number> &Kpv = *elem_subjacobians[p_var][v_var];
  DenseSubMatrix<Number> &Kpw = *elem_subjacobians[p_var][w_var];
  DenseSubVector<Number> &Fp = *elem_subresiduals[p_var];

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      Gradient grad_u = interior_gradient(u_var, qp),
               grad_v = interior_gradient(v_var, qp),
               grad_w = interior_gradient(w_var, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * psi[i][qp] *
                   (grad_u(0) + grad_v(1));
          if (dim == 3)
            Fp(i) += JxW[qp] * psi[i][qp] *
                     (grad_w(2));

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                if (dim == 3)
                  Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2);
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
  const std::vector<std::vector<Real> >& phi_side =
    fe_side_vel->get_phi();

  // The number of local degrees of freedom in u and v
  const unsigned int n_u_dofs = dof_indices_var[u_var].size(); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *elem_subjacobians[w_var][w_var];

  DenseSubVector<Number> &Fu = *elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *elem_subresiduals[w_var];

  // For this example we will use Dirichlet velocity boundary
  // conditions imposed at each timestep via the penalty method.

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const Real penalty = 1.e10;

  unsigned int n_sidepoints = side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      Number u = side_value(u_var, qp),
             v = side_value(v_var, qp),
             w = side_value(w_var, qp);

      // Set u = 1 on the top boundary, (which build_square() has
      // given boundary id 2, and build_cube() has called 5),
      // u = 0 everywhere else
      short int boundary_id =
        this->get_mesh().boundary_info->boundary_id(elem, side);
      assert (boundary_id != BoundaryInfo::invalid_id);
      const short int top_id = (dim==3) ? 5 : 2;

      Real u_value = 0.;

      // For lid-driven cavity, set u=1 on the lid.
      if ((application == 0) && (boundary_id == top_id))
	u_value = 1.;
              
      // Set v, w = 0 everywhere
      const Real v_value = 0.;
      const Real w_value = 0.;

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * penalty *
                   (u - u_value) * phi_side[i][qp];
          Fv(i) += JxW_side[qp] * penalty *
                   (v - v_value) * phi_side[i][qp];
          if (dim == 3)
            Fw(i) += JxW_side[qp] * penalty *
                     (w - w_value) * phi_side[i][qp];

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW_side[qp] * penalty *
                            phi_side[i][qp] * phi_side[j][qp];
                Kvv(i,j) += JxW_side[qp] * penalty *
                            phi_side[i][qp] * phi_side[j][qp];
                if (dim == 3)
                  Kww(i,j) += JxW_side[qp] * penalty *
                              phi_side[i][qp] * phi_side[j][qp];
              }
        }
    }

  // Pin p = 0 at the origin
  if (elem->contains_point(Point(0.,0.)))
    {
      // The pressure penalty value.  \f$ \frac{1}{\epsilon} \f$
      const Real penalty = 1.e9;

      DenseSubMatrix<Number> &Kpp = *elem_subjacobians[p_var][p_var];
      DenseSubVector<Number> &Fp = *elem_subresiduals[p_var];
      const unsigned int n_p_dofs = dof_indices_var[p_var].size(); 

      Point zero(0.);
      Number p = point_value(p_var, zero);
      Number p_value = 0.;

      unsigned int dim = get_mesh().mesh_dimension();
      FEType fe_type = element_fe_var[p_var]->get_fe_type();
      Point p_master = FEInterface::inverse_map(dim, fe_type, elem, zero);

      std::vector<Real> point_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          point_phi[i] = FEInterface::shape(dim, fe_type, elem, i, p_master);
        }

      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += penalty * (p - p_value) * point_phi[i];
          if (request_jacobian)
            for (unsigned int j=0; j != n_p_dofs; j++)
              {
		Kpp(i,j) += penalty * point_phi[i] * point_phi[j];
              }
        }
    }

  return request_jacobian;
}



Point NavierSystem::forcing(const Point& p)
{
  switch (application)
    {
      // lid driven cavity
    case 0:
      {
	// No forcing
	return Point(0.,0.,0.);
      }

      // Homogeneous Dirichlet BCs + sinusoidal forcing
    case 1:
      {
	const unsigned int dim = this->get_mesh().mesh_dimension();

	// This assumes your domain is defined on [0,1]^dim.
	Point f;

	// Counter-Clockwise vortex in x-y plane
	if (dim==2)
	  {
	    f(0) =  std::sin(2.*libMesh::pi*p(1));
	    f(1) = -std::sin(2.*libMesh::pi*p(0));
	  }

	// Counter-Clockwise vortex in x-z plane
	else if (dim==3)
	  {
	    f(0) =  std::sin(2.*libMesh::pi*p(1));
	    f(1) = 0.;
	    f(2) = -std::sin(2.*libMesh::pi*p(0));
	  }

	return f;
      }

    default:
      {
	error();
      }
    }
}
