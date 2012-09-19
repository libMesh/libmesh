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


#include "libmesh/boundary_info.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"


#include "coupled_system.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Boundary conditions for the 3D test case
class BdyFunction : public FunctionBase<Number>
{
public:
  BdyFunction (unsigned int u_var, unsigned int v_var, int sign)
    : _u_var(u_var), _v_var(v_var), _sign(sign)
    { this->_initialized = true; }

  virtual Number operator() (const Point&, const Real = 0)
    { libmesh_not_implemented(); }

  virtual void operator() (const Point& p,
                           const Real,
                           DenseVector<Number>& output)
    {
      output.resize(2);
      output.zero();
      const Real x=p(0), y=p(1);
      output(_u_var) = (_sign)*((y-2) * (y-3));
      output(_v_var) = 0;      
    }

  virtual AutoPtr<FunctionBase<Number> > clone() const
    { return AutoPtr<FunctionBase<Number> > (new BdyFunction(_u_var, _v_var, _sign)); }

private:
  const unsigned int _u_var, _v_var;
  const Real _sign;
};


void CoupledSystem::init_data ()
{  
  // Check the input file for Reynolds number, application type,
  // approximation type
  GetPot infile("coupled_system.in");  
  Peclet = infile("Peclet", 1.);
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
  v_var = this->add_variable ("v", static_cast<Order>(pressure_p+1),
			      fefamily);
  
  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  p_var = this->add_variable ("p", static_cast<Order>(pressure_p),
			      fefamily);

  
  // Tell the system to march velocity forward in time, but 
  // leave p as a constraint only
  this->time_evolving(u_var);
  this->time_evolving(v_var);
  

  // Useful debugging options
  // Set verify_analytic_jacobians to 1e-6 to use
  this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
  this->print_jacobians = infile("print_jacobians", false);
  this->print_element_jacobians = infile("print_element_jacobians", false);

  // Set Dirichlet boundary conditions
  const boundary_id_type left_inlet_id = 0;
  std::set<boundary_id_type> left_inlet_bdy;
  left_inlet_bdy.insert(left_inlet_id);

  const boundary_id_type right_inlet_id = 1;
  std::set<boundary_id_type> right_inlet_bdy;
  right_inlet_bdy.insert(right_inlet_id);

  const boundary_id_type outlets_id = 2;
  std::set<boundary_id_type> outlets_bdy;
  outlets_bdy.insert(outlets_id);

  const boundary_id_type wall_id = 3;
  std::set<boundary_id_type> wall_bdy;
  wall_bdy.insert(wall_id);
  
  std::vector<unsigned int> u_only(1, u_var);
  std::vector<unsigned int> uv(1, u_var);
  uv.push_back(v_var);
      
  // The zero and constant functions
  ZeroFunction<Number> zero;
  ConstFunction<Number> one(1);

  // We need two functions for the inlets, because the signs on the velocities
  // will be different
  int velocity_sign = 1;
  BdyFunction inflow_left(u_var, v_var, -velocity_sign);
  BdyFunction inflow_right(u_var, v_var, velocity_sign);

  // On the walls we will apply the no slip boundary condition, u=0, v=0
  this->get_dof_map().add_dirichlet_boundary
        (DirichletBoundary (wall_bdy, uv, &zero));

  // On the inlet (left)
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (left_inlet_bdy, uv, &inflow_left));
  
   
  // On the inlet (right)
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (right_inlet_bdy, uv, &inflow_right));
  
  //if(species_transport)
  //{
      // Add the Concentration variable "C". They will
      // be approximated using second-order approximation
      C_var = this->add_variable ("C", static_cast<Order>(pressure_p+1),
			      fefamily);
      this->time_evolving(C_var);
      std::vector<unsigned int> C_only(1, C_var);
      this->get_dof_map().add_dirichlet_boundary
	(DirichletBoundary (left_inlet_bdy, C_only, &one));
      this->get_dof_map().add_dirichlet_boundary
	(DirichletBoundary (right_inlet_bdy, C_only, &zero));
      //}
  

  // Do the parent's initialization after variables and boundary constraints are defined
  FEMSystem::init_data();
}



void CoupledSystem::init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system.
  // Note that the concentration and velocity components
  // use the same basis.
  c.element_fe_var[u_var]->get_JxW();
  c.element_fe_var[u_var]->get_phi();
  c.element_fe_var[u_var]->get_dphi();
  c.element_fe_var[u_var]->get_xyz();
  
  c.element_fe_var[p_var]->get_phi();

  c.side_fe_var[u_var]->get_JxW();
  c.side_fe_var[u_var]->get_phi();
  c.side_fe_var[u_var]->get_xyz();
}


bool CoupledSystem::element_time_derivative (bool request_jacobian,
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
  libmesh_assert_equal_to (n_u_dofs, c.dof_indices_var[v_var].size()); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kuv = *c.elem_subjacobians[u_var][v_var];
  DenseSubMatrix<Number> &Kup = *c.elem_subjacobians[u_var][p_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  
  DenseSubMatrix<Number> &Kvu = *c.elem_subjacobians[v_var][u_var];  
  DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];      
  DenseSubMatrix<Number> &Kvp = *c.elem_subjacobians[v_var][p_var];    
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];  

  DenseSubMatrix<Number> &KCu = *c.elem_subjacobians[C_var][u_var];
  DenseSubMatrix<Number> &KCv = *c.elem_subjacobians[C_var][v_var];  
  DenseSubMatrix<Number> &KCC = *c.elem_subjacobians[C_var][C_var];
  DenseSubVector<Number> &FC = *c.elem_subresiduals[C_var];
      
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
	v = c.interior_value(v_var, qp);	
      Gradient grad_u = c.interior_gradient(u_var, qp),
	grad_v = c.interior_gradient(v_var, qp),
	grad_C = c.interior_gradient(C_var, qp);	

      // Definitions for convenience.  It is sometimes simpler to do a
      // dot product if you have the full vector at your disposal.
      NumberVectorValue U     (u,     v);      
      const Number u_x = grad_u(0);
      const Number u_y = grad_u(1);      
      const Number v_x = grad_v(0);
      const Number v_y = grad_v(1);            
      const Number C_x = grad_C(0);
      const Number C_y = grad_C(1);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] *
                   (p*dphi[i][qp](0) -                // pressure term
		    (grad_u*dphi[i][qp]));            // diffusion term                  		    
            
          Fv(i) += JxW[qp] *
                   (p*dphi[i][qp](1) -                // pressure term
		    (grad_v*dphi[i][qp]));            // diffusion term		    		    	    

	  FC(i) += JxW[qp] * ( (U*grad_C)*phi[i][qp] + (1./Peclet)*(grad_C*dphi[i][qp]) ); // Concentration Equation Residual	      	    

          // Note that the Fp block is identically zero unless we are using
          // some kind of artificial compressibility scheme...

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); /* diffusion term  */
                         		  
                  Kvv(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); /* diffusion term  */
		                    
		  KCu(i,j) += JxW[qp]* ( (phi[j][qp]*C_x)*phi[i][qp] );
		  
		  KCv(i,j) += JxW[qp]*( (phi[j][qp]*C_y)*phi[i][qp] );		  	   
		  
		  KCC(i,j) += JxW[qp]*
		    ( (U*dphi[j][qp])*phi[i][qp] + (1./Peclet)*(dphi[j][qp]*dphi[i][qp]) );

		}
	      
	      // Matrix contributions for the up and vp couplings.
	      for (unsigned int j=0; j != n_p_dofs; j++)
		{
		  Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
		  Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
		}
	    }
	}
        
    } // end of the quadrature point qp-loop

      return request_jacobian;    
}


bool CoupledSystem::element_constraint (bool request_jacobian,
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
  DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      Gradient grad_u = c.interior_gradient(u_var, qp),
	grad_v = c.interior_gradient(v_var, qp);	

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * psi[i][qp] *
                   (grad_u(0) + grad_v(1));          

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                  Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);                  
                }
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

void CoupledSystem::postprocess()
{  
  // We need to overload the postprocess function
  // to get the contributions to computed_QoI from each parallel node

  computed_QoI = 0.0;

  FEMSystem::postprocess();

  Parallel::sum(computed_QoI);
}



void CoupledSystem::side_postprocess(DiffContext &context)
{
  
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.side_fe_var[0]->get_JxW();

  const std::vector<Point > &q_point = c.side_fe_var[0]->get_xyz();

  const std::vector<Point> &face_normals = c.side_fe_var[0]->get_normals();

  // Loop over qp's, compute the function at each qp and add
  // to get the QoI

  unsigned int n_qpoints = c.side_qrule->n_points();  
    
  // Get the bdry id for this bdry
  short int bc_id = this->get_mesh().boundary_info->boundary_id (c.elem, c.side);
  if (bc_id==BoundaryInfo::invalid_id)
    {
      std::cout<<"Invalid Boundary id "<<std::endl<<std::endl;
      libmesh_error();
    }

  Real dQoI_0 = 0. ;
  Number u = 0. ;
  Number C = 0. ;
  
  // If side is on the left outlet
  if(bc_id == 2) 
    {             
      //Loop over all the qps on this side 
      for (unsigned int qp=0; qp != n_qpoints; qp++) 
	{ 
	  Real x = q_point[qp](0);

	  if(x < 0.)
	    {
	      // Get u and C at the qp 
	      u = c.side_value(0,qp);
	      C = c.side_value(3,qp);
	  
	      dQoI_0 += JxW[qp] * -u * C;
	    } // end if
	  
  	} // end quadrature loop

    } // end if on bdry
    
  computed_QoI += dQoI_0;
    
}



void CoupledSystem::side_qoi_derivative (DiffContext &context,
                                      const QoISet &qois)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
   
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.side_fe_var[0]->get_JxW();

  // Get velocity basis functions phi
  const std::vector<std::vector<Real> > &phi = c.side_fe_var[0]->get_phi();

  const std::vector<Point > &q_point = c.side_fe_var[0]->get_xyz();
          
  // The number of local degrees of freedom in each variable  
  const unsigned int n_u_dofs = c.dof_indices_var[1].size();  
        
  DenseSubVector<Number> &Qu = *c.elem_qoi_subderivatives[0][0];
  DenseSubVector<Number> &QC = *c.elem_qoi_subderivatives[0][3];
  
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.side_qrule->n_points();  
  
  // Get the bdry id for this bdry
  short int bc_id = this->get_mesh().boundary_info->boundary_id (c.elem, c.side);
  if (bc_id==BoundaryInfo::invalid_id)
    {
      std::cout<<"Invalid Boundary id "<<std::endl<<std::endl;
      libmesh_error();
    }

  Number u = 0. ;
  Number C = 0. ;

  // // If side is on outlet
  if(bc_id == 2)
    {
      // Loop over all the qps on this side
      for (unsigned int qp=0; qp != n_qpoints; qp++)
  	{      
	  Real x = q_point[qp](0);

	  if(x < 0.)
	    {
	      // Get u at the qp 
	      u = c.side_value(0,qp);
	      C = c.side_value(3,qp);
	      
	      // Add the contribution from each basis function
	      for (unsigned int i=0; i != n_u_dofs; i++)
		{
		  Qu(i) += JxW[qp] * -phi[i][qp] * C;
		  QC(i) += JxW[qp] * phi[i][qp] * -u;
		}
	    } // end if
	    	  
  	} // end quadrature loop

    } // end if on outlet    

}



Number CoupledFEMFunctionsx::operator()(const FEMContext& c, const Point& p,
			     const Real time)
{
  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
 
  Real w = 1.;
  
  Number u = c.point_value(0, p);
	
  w = u;		
  
  return w;

}



Number CoupledFEMFunctionsy::operator()(const FEMContext& c, const Point& p,
			     const Real time)
{
  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
 
  Real w = 1.;
  
  Number v = c.point_value(1, p);
	
  w = v;	
	    
  return w;      
}

