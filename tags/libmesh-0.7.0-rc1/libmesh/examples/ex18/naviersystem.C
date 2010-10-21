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
#include "fem_context.h"
#include "mesh.h"
#include "quadrature.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


void NavierSystem::init_data ()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  u_var = this->add_variable ("u", SECOND);
  v_var = this->add_variable ("v", SECOND);

  if (dim == 3)
    w_var = this->add_variable ("w", SECOND);
  else
    w_var = u_var;

  p_var = this->add_variable ("p", FIRST);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Tell the system to march velocity forward in time, but 
  // leave p as a constraint only
  this->time_evolving(u_var);
  this->time_evolving(v_var);
  if (dim == 3)
    this->time_evolving(w_var);

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



void NavierSystem::init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system.
  c.element_fe_var[u_var]->get_JxW();
  c.element_fe_var[u_var]->get_phi();
  c.element_fe_var[u_var]->get_dphi();
  c.element_fe_var[u_var]->get_xyz();
  
  c.element_fe_var[p_var]->get_phi();

  c.side_fe_var[u_var]->get_JxW();
  c.side_fe_var[u_var]->get_phi();
  c.side_fe_var[u_var]->get_xyz();
}


bool NavierSystem::element_time_derivative (bool request_jacobian,
                                            DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = 
    c.element_fe_var[u_var]->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[p_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[u_var]->get_xyz();
 
  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size(); 
  libmesh_assert (n_u_dofs == c.dof_indices_var[v_var].size()); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var];
  DenseSubMatrix<Number> &Kuv = *c.elem_subjacobians[u_var][v_var];
  DenseSubMatrix<Number> &Kuw = *c.elem_subjacobians[u_var][w_var];
  DenseSubMatrix<Number> &Kvu = *c.elem_subjacobians[v_var][u_var];
  DenseSubMatrix<Number> &Kvw = *c.elem_subjacobians[v_var][w_var];
  DenseSubMatrix<Number> &Kwu = *c.elem_subjacobians[w_var][u_var];
  DenseSubMatrix<Number> &Kwv = *c.elem_subjacobians[w_var][v_var];
  DenseSubMatrix<Number> &Kup = *c.elem_subjacobians[u_var][p_var];
  DenseSubMatrix<Number> &Kvp = *c.elem_subjacobians[v_var][p_var];
  DenseSubMatrix<Number> &Kwp = *c.elem_subjacobians[w_var][p_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      Number p = c.interior_value(p_var, qp),
             u = c.interior_value(u_var, qp),
             v = c.interior_value(v_var, qp),
             w = c.interior_value(w_var, qp);
      Gradient grad_u = c.interior_gradient(u_var, qp),
               grad_v = c.interior_gradient(v_var, qp),
               grad_w = c.interior_gradient(w_var, qp);

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

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              // Matrix contributions for the uu and vv couplings.
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
              for (unsigned int j=0; j != n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                  Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
                  if (dim == 3)
                    Kwp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);
                }
            }
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}



bool NavierSystem::element_constraint (bool request_jacobian,
                                       DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[p_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kpu = *c.elem_subjacobians[p_var][u_var];
  DenseSubMatrix<Number> &Kpv = *c.elem_subjacobians[p_var][v_var];
  DenseSubMatrix<Number> &Kpw = *c.elem_subjacobians[p_var][w_var];
  DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      Gradient grad_u = c.interior_gradient(u_var, qp),
               grad_v = c.interior_gradient(v_var, qp),
               grad_w = c.interior_gradient(w_var, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * psi[i][qp] *
                   (grad_u(0) + grad_v(1));
          if (dim == 3)
            Fp(i) += JxW[qp] * psi[i][qp] *
                     (grad_w(2));

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                  Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                  if (dim == 3)
                    Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2);
                }
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

      

bool NavierSystem::side_constraint (bool request_jacobian,
                                    DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = c.side_fe_var[u_var]->get_JxW();

  // The velocity shape functions at side quadrature points.
  const std::vector<std::vector<Real> >& phi_side =
    c.side_fe_var[u_var]->get_phi();

  // Physical location of the quadrature points on the side
  const std::vector<Point>& qpoint = c.side_fe_var[u_var]->get_xyz();

  // The number of local degrees of freedom in u and v
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size(); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var];

  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];

  // For this example we will use Dirichlet velocity boundary
  // conditions imposed at each timestep via the penalty method.

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const Real penalty = 1.e10;

  unsigned int n_sidepoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      Number u = c.side_value(u_var, qp),
             v = c.side_value(v_var, qp),
             w = c.side_value(w_var, qp);

      // Set u = 1 on the top boundary, (which build_square() has
      // given boundary id 2, and build_cube() has called 5),
      // u = 0 everywhere else
      short int boundary_id =
        this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
      libmesh_assert (boundary_id != BoundaryInfo::invalid_id);
      const short int top_id = (dim==3) ? 5 : 2;

      // Boundary data coming from true solution
      Real
	u_value = 0.,
	v_value = 0.,
	w_value = 0.;

      // For lid-driven cavity, set u=1 on the lid.
      if ((application == 0) && (boundary_id == top_id))
	u_value = 1.;

      if (application == 2)
	{
	  const Real x=qpoint[qp](0);
	  const Real y=qpoint[qp](1);
	  const Real z=qpoint[qp](2);
	  u_value=(Reynolds+1)*(y*y + z*z);
	  v_value=(Reynolds+1)*(x*x + z*z);
	  w_value=(Reynolds+1)*(x*x + y*y);
	}
      
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * penalty *
                   (u - u_value) * phi_side[i][qp];
          Fv(i) += JxW_side[qp] * penalty *
                   (v - v_value) * phi_side[i][qp];
          if (dim == 3)
            Fw(i) += JxW_side[qp] * penalty *
                     (w - w_value) * phi_side[i][qp];

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

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
    }

  // Pin p = 0 at the origin
  if (c.elem->contains_point(Point(0.,0.)))
    {
      // The pressure penalty value.  \f$ \frac{1}{\epsilon} \f$
      const Real penalty = 1.e9;

      DenseSubMatrix<Number> &Kpp = *c.elem_subjacobians[p_var][p_var];
      DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size(); 

      Point zero(0.);
      Number p = c.point_value(p_var, zero);
      Number p_value = 0.;

      unsigned int dim = get_mesh().mesh_dimension();
      FEType fe_type = c.element_fe_var[p_var]->get_fe_type();
      Point p_master = FEInterface::inverse_map(dim, fe_type, c.elem, zero);

      std::vector<Real> point_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          point_phi[i] = FEInterface::shape(dim, fe_type, c.elem, i, p_master);
        }

      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += penalty * (p - p_value) * point_phi[i];
          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_p_dofs; j++)
		Kpp(i,j) += penalty * point_phi[i] * point_phi[j];
            }
        }
    }

  return request_jacobian;
}


// We override the default mass_residual function,
// because in the non-dimensionalized Navier-Stokes equations
// the time derivative of velocity has a Reynolds number coefficient.
// Alternatively we could divide the whole equation by
// Reynolds number (or choose a more complicated non-dimensionalization
// of time), but this gives us an opportunity to demonstrate overriding
// FEMSystem::mass_residual()
bool NavierSystem::mass_residual (bool request_jacobian,
                                  DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = 
    c.element_fe_var[u_var]->get_phi();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var];

  // The number of local degrees of freedom in velocity
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      Number u = c.interior_value(u_var, qp),
             v = c.interior_value(v_var, qp),
             w = c.interior_value(w_var, qp);

      // We pull as many calculations as possible outside of loops
      Number JxWxRe   = JxW[qp] * Reynolds;
      Number JxWxRexU = JxWxRe * u;
      Number JxWxRexV = JxWxRe * v;
      Number JxWxRexW = JxWxRe * w;

      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
          Fu(i) += JxWxRexU * phi[i][qp];
          Fv(i) += JxWxRexV * phi[i][qp];
          if (dim == 3)
            Fw(i) += JxWxRexW * phi[i][qp];

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              Number JxWxRexPhiI = JxWxRe * phi[i][qp];
              Number JxWxRexPhiII = JxWxRexPhiI * phi[i][qp];
              Kuu(i,i) += JxWxRexPhiII;
              Kvv(i,i) += JxWxRexPhiII;
              if (dim == 3)
                Kww(i,i) += JxWxRexPhiII;

              // The mass matrix is symmetric, so we calculate
              // one triangle and add it to both upper and lower
              // triangles
              for (unsigned int j = i+1; j != n_u_dofs; ++j)
                {
                  Number Kij = JxWxRexPhiI * phi[j][qp];
                  Kuu(i,j) += Kij;
                  Kuu(j,i) += Kij;
                  Kvv(i,j) += Kij;
                  Kvv(j,i) += Kij;
                  if (dim == 3)
                    {
                      Kww(i,j) += Kij;
                      Kww(j,i) += Kij;
                    }
                }
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

      // 3D test case with quadratic velocity and linear pressure field
    case 2:
      {
	// This problem doesn't make sense in 1D...
	libmesh_assert(1 != this->get_mesh().mesh_dimension());
	const Real x=p(0), y=p(1), z=p(2);
	const Real
	  u=(Reynolds+1)*(y*y + z*z),
	  v=(Reynolds+1)*(x*x + z*z),
	  w=(Reynolds+1)*(x*x + y*y);

	if (this->get_mesh().mesh_dimension() == 2)
	  return 2.*Reynolds*(Reynolds+1.)*Point(v*y,
						 u*x);
	else
	  return 2.*Reynolds*(Reynolds+1.)*Point(v*y + w*z,
						 u*x + w*z,
						 u*x + v*y);
      }

    default:
      {
	libmesh_error();
      }
    }
}
