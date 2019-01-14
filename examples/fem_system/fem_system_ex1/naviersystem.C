// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/getpot.h"

#include "naviersystem.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"
#include "libmesh/elem.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Boundary conditions for the 3D test case
class BdyFunction : public FunctionBase<Number>
{
public:
  BdyFunction (unsigned int u_var,
               unsigned int v_var,
               unsigned int w_var,
               Real Reynolds)
    : _u_var(u_var), _v_var(v_var), _w_var(w_var), _Re(Reynolds)
  { this->_initialized = true; }

  virtual Number operator() (const Point &, const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output)
  {
    output.zero();
    const Real x=p(0), y=p(1), z=p(2);
    output(_u_var) = (_Re+1)*(y*y + z*z);
    output(_v_var) = (_Re+1)*(x*x + z*z);
    output(_w_var) = (_Re+1)*(x*x + y*y);
  }

  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  { return libmesh_make_unique<BdyFunction>(_u_var, _v_var, _w_var, _Re); }

private:
  const unsigned int _u_var, _v_var, _w_var;
  const Real _Re;
};


void NavierSystem::init_data ()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Check the input file for Reynolds number, application type,
  // approximation type
  GetPot infile("navier.in");
  Reynolds = infile("Reynolds", 1.);
  application = infile("application", 0);
  unsigned int pressure_p = infile("pressure_p", 1);
  std::string fe_family = infile("fe_family", std::string("LAGRANGE"));

  // LBB needs better-than-quadratic velocities for better-than-linear
  // pressures, and libMesh needs non-Lagrange elements for
  // better-than-quadratic velocities.
  libmesh_assert((pressure_p == 1) || (fe_family != "LAGRANGE"));

  FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);

  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  u_var = this->add_variable ("u", static_cast<Order>(pressure_p+1),
                              fefamily);
  v_var = this->add_variable ("v",
                              static_cast<Order>(pressure_p+1),
                              fefamily);

  if (dim == 3)
    w_var = this->add_variable ("w",
                                static_cast<Order>(pressure_p+1),
                                fefamily);
  else
    w_var = u_var;

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  p_var = this->add_variable ("p",
                              static_cast<Order>(pressure_p),
                              fefamily);

  // Tell the system to march velocity forward in time, but
  // leave p as a constraint only
  this->time_evolving(u_var, 1);
  this->time_evolving(v_var, 1);
  if (dim == 3)
    this->time_evolving(w_var, 1);

  // Useful debugging options
  // Set verify_analytic_jacobians to 1e-6 to use
  this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
  this->print_jacobians = infile("print_jacobians", false);
  this->print_element_jacobians = infile("print_element_jacobians", false);

  // Set Dirichlet boundary conditions
  const boundary_id_type top_id = (dim==3) ? 5 : 2;

  std::set<boundary_id_type> top_bdys;
  top_bdys.insert(top_id);

  const boundary_id_type all_ids[6] = {0, 1, 2, 3, 4, 5};
  std::set<boundary_id_type> all_bdys(all_ids, all_ids+(dim*2));

  std::set<boundary_id_type> nontop_bdys = all_bdys;
  nontop_bdys.erase(top_id);

  std::vector<unsigned int> u_only(1, u_var);
  std::vector<unsigned int> vw(1, v_var), uvw(1, u_var);
  uvw.push_back(v_var);
  if (dim == 3)
    {
      vw.push_back(w_var);
      uvw.push_back(w_var);
    }

  ZeroFunction<Number> zero;
  ConstFunction<Number> one(1);
  // For lid-driven cavity, set u=1,v=w=0 on the lid and u=v=w=0 elsewhere
  if (application == 0)
    {
      this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (top_bdys, u_only, &one));
      this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (top_bdys, vw, &zero));
      this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (nontop_bdys, uvw, &zero));
    }
  // For forcing with zero wall velocity, set homogeneous Dirichlet BCs
  else if (application == 1)
    {
      this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (all_bdys, uvw, &zero));
    }
  // For 3D test case with quadratic velocity field, set that field on walls
  else if (application == 2)
    {
      BdyFunction bdy(u_var, v_var, w_var, Reynolds);
      this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (all_bdys, uvw, &bdy));
    }

  // Do the parent's initialization after variables and boundary constraints are defined
  FEMSystem::init_data();
}



void NavierSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * u_elem_fe;
  FEBase * p_elem_fe;
  FEBase * u_side_fe;

  c.get_element_fe(u_var, u_elem_fe);
  c.get_element_fe(p_var, p_elem_fe);
  c.get_side_fe(u_var, u_side_fe);

  // We should prerequest all the data
  // we will need to build the linear system.
  u_elem_fe->get_JxW();
  u_elem_fe->get_phi();
  u_elem_fe->get_dphi();
  u_elem_fe->get_xyz();

  p_elem_fe->get_phi();

  u_side_fe->get_JxW();
  u_side_fe->get_phi();
  u_side_fe->get_xyz();
}


bool NavierSystem::element_time_derivative (bool request_jacobian,
                                            DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * u_elem_fe;
  FEBase * p_elem_fe;

  c.get_element_fe(u_var, u_elem_fe);
  c.get_element_fe(p_var, p_elem_fe);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = u_elem_fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real>> & phi = u_elem_fe->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = u_elem_fe->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real>> & psi = p_elem_fe->get_phi();

  // Physical location of the quadrature points
  const std::vector<Point> & qpoint = u_elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = c.n_dof_indices(p_var);
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);
  libmesh_assert_equal_to (n_u_dofs, c.n_dof_indices(v_var));

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(u_var, u_var);
  DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);
  DenseSubMatrix<Number> & Kww = c.get_elem_jacobian(w_var, w_var);
  DenseSubMatrix<Number> & Kuv = c.get_elem_jacobian(u_var, v_var);
  DenseSubMatrix<Number> & Kuw = c.get_elem_jacobian(u_var, w_var);
  DenseSubMatrix<Number> & Kvu = c.get_elem_jacobian(v_var, u_var);
  DenseSubMatrix<Number> & Kvw = c.get_elem_jacobian(v_var, w_var);
  DenseSubMatrix<Number> & Kwu = c.get_elem_jacobian(w_var, u_var);
  DenseSubMatrix<Number> & Kwv = c.get_elem_jacobian(w_var, v_var);
  DenseSubMatrix<Number> & Kup = c.get_elem_jacobian(u_var, p_var);
  DenseSubMatrix<Number> & Kvp = c.get_elem_jacobian(v_var, p_var);
  DenseSubMatrix<Number> & Kwp = c.get_elem_jacobian(w_var, p_var);
  DenseSubVector<Number> & Fu = c.get_elem_residual(u_var);
  DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);
  DenseSubVector<Number> & Fw = c.get_elem_residual(w_var);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Variables to store solutions & its gradient at old Newton iterate
  Number p, u, v, w;
  Gradient grad_u, grad_v, grad_w;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      c.interior_value(p_var, qp, p);
      c.interior_value(u_var, qp, u);
      c.interior_value(v_var, qp, v);
      c.interior_value(w_var, qp, w);
      c.interior_gradient(u_var, qp, grad_u);
      c.interior_gradient(v_var, qp, grad_v);
      c.interior_gradient(w_var, qp, grad_w);

      // Definitions for convenience.  It is sometimes simpler to do a
      // dot product if you have the full vector at your disposal.
      NumberVectorValue U (u, v);
      if (dim == 3)
        U(2) = w;
      const Number u_x = grad_u(0);
      const Number u_y = grad_u(1);
      const Number u_z = (dim == 3) ? grad_u(2) : 0;
      const Number v_x = grad_v(0);
      const Number v_y = grad_v(1);
      const Number v_z = (dim == 3) ? grad_v(2) : 0;
      const Number w_x = (dim == 3) ? grad_w(0) : 0;
      const Number w_y = (dim == 3) ? grad_w(1) : 0;
      const Number w_z = (dim == 3) ? grad_w(2) : 0;

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
              libmesh_assert_equal_to (c.elem_solution_derivative, 1.0);

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
                                       DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * u_elem_fe;
  FEBase * p_elem_fe;

  c.get_element_fe(u_var, u_elem_fe);
  c.get_element_fe(p_var, p_elem_fe);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> & JxW = u_elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = u_elem_fe->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real>> & psi = p_elem_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);
  const unsigned int n_p_dofs = c.n_dof_indices(p_var);

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> & Kpu = c.get_elem_jacobian(p_var, u_var);
  DenseSubMatrix<Number> & Kpv = c.get_elem_jacobian(p_var, v_var);
  DenseSubMatrix<Number> & Kpw = c.get_elem_jacobian(p_var, w_var);
  DenseSubVector<Number> & Fp = c.get_elem_residual(p_var);

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Variables to store solutions & its gradient at old Newton iterate
  Gradient grad_u, grad_v, grad_w;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      c.interior_gradient(u_var, qp, grad_u),
        c.interior_gradient(v_var, qp, grad_v),
        c.interior_gradient(w_var, qp, grad_w);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * psi[i][qp] *
            (grad_u(0) + grad_v(1));
          if (dim == 3)
            Fp(i) += JxW[qp] * psi[i][qp] *
              (grad_w(2));

          if (request_jacobian && c.get_elem_solution_derivative())
            {
              libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);

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
                                    DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * p_elem_fe;

  c.get_element_fe(p_var, p_elem_fe);

  // Pin p = 0 at the origin
  const Point zero(0., 0.);

  if (c.get_elem().contains_point(zero))
    {
      // The pressure penalty value.  \f$ \frac{1}{\epsilon} \f$
      const Real penalty = 1.e9;

      DenseSubMatrix<Number> & Kpp = c.get_elem_jacobian(p_var, p_var);
      DenseSubVector<Number> & Fp = c.get_elem_residual(p_var);
      const unsigned int n_p_dofs = c.n_dof_indices(p_var);

      Number p;
      c.point_value(p_var, zero, p);

      Number p_value = 0.;

      unsigned int dim = get_mesh().mesh_dimension();
      FEType fe_type = p_elem_fe->get_fe_type();
      Point p_master = FEInterface::inverse_map(dim, fe_type, &c.get_elem(), zero);

      std::vector<Real> point_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          point_phi[i] = FEInterface::shape(dim, fe_type, &c.get_elem(), i, p_master);
        }

      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += penalty * (p - p_value) * point_phi[i];
          if (request_jacobian && c.get_elem_solution_derivative())
            {
              libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);

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
                                  DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * u_elem_fe;

  c.get_element_fe(u_var, u_elem_fe);

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> & JxW = u_elem_fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real>> & phi = u_elem_fe->get_phi();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Number> & Fu = c.get_elem_residual(u_var);
  DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);
  DenseSubVector<Number> & Fw = c.get_elem_residual(w_var);
  DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(u_var, u_var);
  DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);
  DenseSubMatrix<Number> & Kww = c.get_elem_jacobian(w_var, w_var);

  // The number of local degrees of freedom in velocity
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Variables to store time derivatives at old Newton iterate
  Number u_dot, v_dot, w_dot;

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // Compute time derivatives
      c.interior_rate(u_var, qp, u_dot),
        c.interior_rate(v_var, qp, v_dot),
        c.interior_rate(w_var, qp, w_dot);

      // We pull as many calculations as possible outside of loops
      Number JxWxRe   = JxW[qp] * Reynolds;
      Number JxWxRexU = JxWxRe * u_dot;
      Number JxWxRexV = JxWxRe * v_dot;
      Number JxWxRexW = JxWxRe * w_dot;

      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
          Fu(i) -= JxWxRexU * phi[i][qp];
          Fv(i) -= JxWxRexV * phi[i][qp];
          if (dim == 3)
            Fw(i) -= JxWxRexW * phi[i][qp];

          if (request_jacobian && c.get_elem_solution_derivative())
            {
              libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);

              Number JxWxRexPhiI = JxWxRe * phi[i][qp];
              Number JxWxRexPhiII = -JxWxRexPhiI * phi[i][qp];
              Kuu(i,i) += JxWxRexPhiII;
              Kvv(i,i) += JxWxRexPhiII;
              if (dim == 3)
                Kww(i,i) += JxWxRexPhiII;

              // The mass matrix is symmetric, so we calculate
              // one triangle and add it to both upper and lower
              // triangles
              for (unsigned int j = i+1; j != n_u_dofs; ++j)
                {
                  Number Kij = -JxWxRexPhiI * phi[j][qp];
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



void NavierSystem::postprocess()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  Point p(1./3., 1./3.);
  Number u = point_value(u_var, p),
    v = point_value(v_var, p),
    w = (dim == 3) ? point_value(w_var, p) : 0;

  libMesh::out << "u(1/3,1/3) = ("
               << u << ", "
               << v << ", "
               << w << ")"
               << std::endl;
}




Point NavierSystem::forcing(const Point & p)
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
        libmesh_assert_not_equal_to (1, this->get_mesh().mesh_dimension());
        const Real x=p(0), y=p(1), z=p(2);
        const Real
          u=(Reynolds+1)*(y*y + z*z),
          v=(Reynolds+1)*(x*x + z*z),
          w=(Reynolds+1)*(x*x + y*y);

        if (this->get_mesh().mesh_dimension() == 2)
          return 2*Reynolds*(Reynolds+1)*Point(v*y,
                                               u*x);
        else
          return 2*Reynolds*(Reynolds+1)*Point(v*y + w*z,
                                               u*x + w*z,
                                               u*x + v*y);
      }

    default:
      libmesh_error_msg("Invalid application id = " << application);
    }
}
