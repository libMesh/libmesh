/* The libMesh Finite Element Library. */
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

#include "libmesh/getpot.h"

#include "curl_curl_system.h"

#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

CurlCurlSystem::CurlCurlSystem( EquationSystems& es,
			      const std::string& name_in,
			      const unsigned int number_in)
  : FEMSystem(es, name_in, number_in)
{
  return;
}

void CurlCurlSystem::init_data ()
{
  // Check the input file for Reynolds number, application type,
  // approximation type
  GetPot infile("vector_fe_ex4.in");

  // Add the solution variable
  u_var = this->add_variable ("u", FIRST, NEDELEC_ONE);

  this->time_evolving(u_var);

  // Useful debugging options
  // Set verify_analytic_jacobians to 1e-6 to use
  this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
  this->print_jacobians = infile("print_jacobians", false);
  this->print_element_jacobians = infile("print_element_jacobians", false);

  this->extra_quadrature_order = infile("extra_quadrature_order", 0 );

  this->init_dirichlet_bcs();

  // Do the parent's initialization after variables and boundary constraints are defined
  FEMSystem::init_data();
}

void CurlCurlSystem::init_dirichlet_bcs()
{
  const boundary_id_type all_ids[6] = {0,1,2,3,4,5};
  std::set<boundary_id_type> boundary_ids( all_ids, all_ids+6 );

  std::vector<unsigned int> vars;
  vars.push_back( u_var );

  ZeroFunction<Real> func;

  // HCurl DirichletBoundary not supported yet
  //this->get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( boundary_ids, vars, &func ) );
  return;
}

void CurlCurlSystem::init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Get finite element object
  FEGenericBase<RealGradient>* fe;
  c.get_element_fe<RealGradient>( u_var, fe );

  // We should prerequest all the data
  // we will need to build the linear system.
  fe->get_JxW();
  fe->get_phi();
  fe->get_curl_phi();
  fe->get_xyz();

  // Get finite element object
  FEGenericBase<RealGradient>* side_fe;
  c.get_side_fe<RealGradient>( u_var, side_fe );

  side_fe->get_JxW();
  side_fe->get_phi();

  return;
}


bool CurlCurlSystem::element_time_derivative (bool request_jacobian,
					      DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Get finite element object
  FEGenericBase<RealGradient>* fe = NULL;
  c.get_element_fe<RealGradient>( u_var, fe );

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<RealGradient> >& phi = fe->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& curl_phi = fe->get_curl_phi();

  const std::vector<Point>& qpoint = fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

  DenseSubMatrix<Number> &Kuu = c.get_elem_jacobian(u_var,u_var);

  DenseSubVector<Number> &Fu = c.get_elem_residual(u_var);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  const unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Loop over quadrature points
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Gradient u;
      Gradient curl_u;

      c.interior_value( u_var, qp, u );

      c.interior_curl( u_var, qp, curl_u );

      // Value of the forcing function at this quadrature point
      RealGradient f = this->forcing(qpoint[qp]);

      // First, an i-loop over the degrees of freedom.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += ( curl_u*curl_phi[i][qp] + u*phi[i][qp] - f*phi[i][qp] )*JxW[qp];

          if (request_jacobian)
            {
              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kuu(i,j) += ( curl_phi[j][qp]*curl_phi[i][qp] +
			        phi[j][qp]*phi[i][qp] )*JxW[qp];

		}
	    }

        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}


bool CurlCurlSystem::side_time_derivative (bool request_jacobian,
					   DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Get finite element object
  FEGenericBase<RealGradient>* side_fe = NULL;
  c.get_side_fe<RealGradient>( u_var, side_fe );

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = side_fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<RealGradient> >& phi = side_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

  const std::vector<Point>& normals = side_fe->get_normals();

  const std::vector<Point>& qpoint = side_fe->get_xyz();

  // The penalty value.  \frac{1}{\epsilon}
  // in the discussion above.
  const Real penalty = 1.e10;

  DenseSubMatrix<Number> &Kuu = c.get_elem_jacobian(u_var,u_var);
  DenseSubVector<Number> &Fu = c.get_elem_residual(u_var);

  const unsigned int n_qpoints = c.get_side_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Gradient u;
      c.side_value( u_var, qp, u );

      RealGradient N( normals[qp](0), normals[qp](1), normals[qp](2) );

      Gradient u_exact = this->exact_solution( qpoint[qp] ) ;

      Gradient Ncu = (u - u_exact).cross(N);

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  Fu(i) += penalty*Ncu*(phi[i][qp].cross(N))*JxW[qp];

          if (request_jacobian)
            {
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kuu(i,j) += penalty*(phi[j][qp].cross(N))*(phi[i][qp].cross(N))*JxW[qp];
            }
        }

    }

  return request_jacobian;
}

RealGradient CurlCurlSystem::forcing( const Point& p )
{
  Real x = p(0); Real y = p(1); Real z = p(2);

  return soln.forcing(x,y,z);
}

RealGradient CurlCurlSystem::exact_solution( const Point& p )
{
  Real x = p(0); Real y = p(1); Real z = p(2);

  return soln(x,y,z);
}
