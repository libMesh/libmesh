/* This is where we define the assembly of the Laplace system */

// General libMesh includes
#include "getpot.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "string_to_enum.h"
#include "equation_systems.h"
#include "parallel.h"
#include "fem_context.h"
#include "femparameters.h"

// Local includes
#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void LaplaceSystem::init_data ()
{
  this->add_variable ("T", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));

  GetPot infile("l-shaped.in");
  parameters.push_back(infile("alpha_1", 1.0));
  parameters.push_back(infile("alpha_2", 1.0));
  
  exact_QoI[0] = infile("QoI_0", 0.0);
  exact_QoI[1] = infile("QoI_1", 0.0);
  
  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  this->time_evolving(0);
}

void LaplaceSystem::init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  c.element_fe_var[0]->get_JxW();
  c.element_fe_var[0]->get_phi();
  c.element_fe_var[0]->get_dphi();

  c.side_fe_var[0]->get_JxW();
  c.side_fe_var[0]->get_phi();
  c.side_fe_var[0]->get_dphi();
}

#define optassert(X) {if (!(X)) libmesh_error();}

// Assemble the element contributions to the stiffness matrix
bool LaplaceSystem::element_time_derivative (bool request_jacobian,
					  DiffContext &context)
{
  // Are the jacobians specified analytically ?
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[0]->get_JxW();

  // Element basis function gradients
  const std::vector<std::vector<RealGradient> > &dphi = c.element_fe_var[0]->get_dphi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[0].size(); 

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &K = *c.elem_subjacobians[0][0];
  DenseSubVector<Number> &F = *c.elem_subresiduals[0];

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {     	
      // Compute the solution gradient at the Newton iterate
      Gradient grad_T = c.interior_gradient(0, qp);

      // The residual contribution from this element
      for (unsigned int i=0; i != n_T_dofs; i++)
        F(i) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( grad_T * dphi[i][qp] ) ;
      if (compute_jacobian)
        for (unsigned int i=0; i != n_T_dofs; i++)
          for (unsigned int j=0; j != n_T_dofs; ++j)
	    // The analytic jacobian
            K(i,j) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( dphi[i][qp] * dphi[j][qp] );
    } // end of the quadrature point qp-loop
  
  return compute_jacobian;
}

// Set Dirichlet bcs, side contributions to global stiffness matrix
bool LaplaceSystem::side_constraint (bool request_jacobian,
				  DiffContext &context)
{
  // Are the jacobians specified analytically ? 
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.side_fe_var[0]->get_JxW();

  // Side basis functions
  const std::vector<std::vector<Real> > &phi = c.side_fe_var[0]->get_phi();

  // Side Quadrature points
  const std::vector<Point > &qside_point = c.side_fe_var[0]->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[0].size(); 

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &K = *c.elem_subjacobians[0][0];
  DenseSubVector<Number> &F = *c.elem_subresiduals[0];
      
  unsigned int n_qpoints = c.side_qrule->n_points();

  Real penalty = 1.e10;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      Number T = c.side_value(0, qp);
            
      // We get the Dirichlet bcs from the exact solution
      Number u_dirichlet = exact_solution (qside_point[qp]);

      // The residual from the boundary terms, penalize non-zero temperature
      for (unsigned int i=0; i != n_T_dofs; i++)
	F(i) += JxW[qp] * penalty * ( T - u_dirichlet) * phi[i][qp];
      if (compute_jacobian)
	for (unsigned int i=0; i != n_T_dofs; i++)
	  for (unsigned int j=0; j != n_T_dofs; ++j)
	    // The analytic jacobian
	    K(i,j) += JxW[qp] * penalty * phi[i][qp] * phi[j][qp];    
	
    } // end of the quadrature point qp-loop
  
  return compute_jacobian;
}

// Override the default DiffSystem postprocess function to compute the 
// approximations to the QoIs 
void LaplaceSystem::postprocess()
{  
  // Reset the array holding the computed QoIs
  computed_QoI[0] = 0.0;
  computed_QoI[1] = 0.0;

  FEMSystem::postprocess();

  Parallel::sum(computed_QoI[0]);

  Parallel::sum(computed_QoI[1]);
     
}

// The exact solution to the singular problem,
// u_exact = r^(2/3)*sin(2*theta/3). We use this to set the Dirichlet boundary conditions
Number LaplaceSystem::exact_solution(const Point& p)// xyz location   
{
  const Real x1 = p(0);
  const Real x2 = p(1);

  Real theta = atan2(x2,x1);

  // Make sure 0 <= theta <= 2*pi
  if (theta < 0)
    theta += 2. * libMesh::pi;
    
  return (parameters[0] + (2.*parameters[1])) * pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
}
