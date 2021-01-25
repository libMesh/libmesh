// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Example includes
#include "naviersystem.h"

// libMesh includes
#include "libmesh/getpot.h"
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

// C++ includes
#include <array>
#include <functional> // std::reference_wrapper

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
#ifdef LIBMESH_ENABLE_DIRICHLET
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
#endif // LIBMESH_ENABLE_DIRICHLET

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

  // And tell the context what data we *don't* need
  FEBase * p_side_fe;
  c.get_side_fe(p_var, p_side_fe);
  p_side_fe->get_nothing();
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

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Get DenseSubMatrix references for velocity-velocity coupling
  // Kuu, Kuv, Kuw
  // Kvu, Kvv, Kvw
  // Kwu, Kwv, Kww
  std::reference_wrapper<DenseSubMatrix<Number>> K[3][3] =
    {
      {c.get_elem_jacobian(u_var, u_var), c.get_elem_jacobian(u_var, v_var), c.get_elem_jacobian(u_var, w_var)},
      {c.get_elem_jacobian(v_var, u_var), c.get_elem_jacobian(v_var, v_var), c.get_elem_jacobian(v_var, w_var)},
      {c.get_elem_jacobian(w_var, u_var), c.get_elem_jacobian(w_var, v_var), c.get_elem_jacobian(w_var, w_var)}
    };

  // Get DenseSubMatrix references for velocity-pressure coupling
  // Kup
  // Kvp
  // Kwp
  std::reference_wrapper<DenseSubMatrix<Number>> B[3] =
    {
      c.get_elem_jacobian(u_var, p_var),
      c.get_elem_jacobian(v_var, p_var),
      c.get_elem_jacobian(w_var, p_var)
    };

  // Get rhs DenseSubVector references
  // Fu
  // Fv
  // Fw
  std::reference_wrapper<DenseSubVector<Number>> F[3] =
    {
      c.get_elem_residual(u_var),
      c.get_elem_residual(v_var),
      c.get_elem_residual(w_var)
    };

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Pressure at current Newton iterate
  Number p;

  // Velocity at current Newton iterate
  NumberVectorValue U;

  // Store (grad_u, grad_v, grad_w) at current Newton iterate
  Gradient grad_uvw[3];

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the current Newton iterate
      c.interior_value(p_var, qp, p);
      c.interior_value(u_var, qp, U(0));
      c.interior_value(v_var, qp, U(1));
      c.interior_value(w_var, qp, U(2));
      c.interior_gradient(u_var, qp, grad_uvw[0]);
      c.interior_gradient(v_var, qp, grad_uvw[1]);
      c.interior_gradient(w_var, qp, grad_uvw[2]);

      // Value of the forcing function at this quadrature point
      Point f = this->forcing(qpoint[qp]);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          for (unsigned int d = 0; d < dim; ++d)
            F[d](i) += JxW[qp] *
              (-Reynolds*(U*grad_uvw[d])*phi[i][qp] + // convection term
               p*dphi[i][qp](d) -                     // pressure term
               (grad_uvw[d]*dphi[i][qp]) +            // diffusion term
               f(d)*phi[i][qp]);                      // forcing function

          // Note that the Fp block is identically zero unless we are using
          // some kind of artificial compressibility scheme...
          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert_equal_to (c.elem_solution_derivative, 1.0);

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j != n_u_dofs; j++)
                for (unsigned int k = 0; k < dim; ++k)
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      // "Diagonal" contributions
                      if (k == l)
                        K[k][k](i,j) += JxW[qp] *
                          (-Reynolds*(U*dphi[j][qp])*phi[i][qp] - // convection term
                           (dphi[i][qp]*dphi[j][qp]));            // diffusion term

                      // Newton terms
                      K[k][l](i,j) += JxW[qp] * -Reynolds*grad_uvw[k](l)*phi[i][qp]*phi[j][qp];
                    }

              // Matrix contributions for the velocity-pressure couplings.
              for (unsigned int j=0; j != n_p_dofs; j++)
                for (unsigned int k = 0; k < dim; ++k)
                  B[k](i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](k);
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

  // Get DenseSubMatrix references for velocity-pressure coupling
  // Kpu Kpv Kpw
  std::reference_wrapper<DenseSubMatrix<Number>> B[3] =
    {
      c.get_elem_jacobian(p_var, u_var),
      c.get_elem_jacobian(p_var, v_var),
      c.get_elem_jacobian(p_var, w_var)
    };

  // Get reference to pressure equation elem residual
  DenseSubVector<Number> & Fp = c.get_elem_residual(p_var);

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Store (grad_u, grad_v, grad_w) at current Newton iterate
  Gradient grad_uvw[3];

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      c.interior_gradient(u_var, qp, grad_uvw[0]);
      c.interior_gradient(v_var, qp, grad_uvw[1]);
      c.interior_gradient(w_var, qp, grad_uvw[2]);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          for (unsigned int k = 0; k < dim; ++k)
            Fp(i) += JxW[qp] * psi[i][qp] * grad_uvw[k](k);

          if (request_jacobian && c.get_elem_solution_derivative())
            {
              libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                for (unsigned int k = 0; k < dim; ++k)
                  B[k](i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](k);
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
      Point p_master = FEMap::inverse_map(dim, &c.get_elem(), zero);

      std::vector<Real> point_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          point_phi[i] = FEInterface::shape(fe_type, &c.get_elem(), i, p_master);
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

  // Get rhs DenseSubVector references
  std::reference_wrapper<DenseSubVector<Number>> F[3] =
    {
      c.get_elem_residual(u_var),
      c.get_elem_residual(v_var),
      c.get_elem_residual(w_var)
    };

  // Get references to stiffness matrix diagonal blocks
  std::reference_wrapper<DenseSubMatrix<Number>> Kdiag[3] =
    {
      c.get_elem_jacobian(u_var, u_var),
      c.get_elem_jacobian(v_var, v_var),
      c.get_elem_jacobian(w_var, w_var)
    };

  // The number of local degrees of freedom in velocity
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Stores JxW * Reynolds * du/dt at current Newton iterate
  std::array<Number, 3> accel;

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // Compute time derivatives
      c.interior_rate(u_var, qp, accel[0]);
      c.interior_rate(v_var, qp, accel[1]);
      c.interior_rate(w_var, qp, accel[2]);

      // We pull as many calculations as possible outside of loops
      Number JxWxRe = JxW[qp] * Reynolds;
      for (unsigned int k = 0; k < dim; ++k)
        accel[k] *= JxWxRe;

      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            F[k](i) -= accel[k] * phi[i][qp];

          if (request_jacobian && c.get_elem_solution_derivative())
            {
              libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);

              Number JxWxRexPhiI = JxWxRe * phi[i][qp];
              Number JxWxRexPhiII = -JxWxRexPhiI * phi[i][qp];
              for (unsigned int k = 0; k < dim; ++k)
                Kdiag[k](i,i) += JxWxRexPhiII;

              // The mass matrix is symmetric, so we calculate
              // one triangle and add it to both upper and lower
              // triangles
              for (unsigned int j = i+1; j != n_u_dofs; ++j)
                {
                  Number Kij = -JxWxRexPhiI * phi[j][qp];

                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      Kdiag[k](i,j) += Kij;
                      Kdiag[k](j,i) += Kij;
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
