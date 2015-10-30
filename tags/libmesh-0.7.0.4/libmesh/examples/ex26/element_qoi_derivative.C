// General libMesh includes
#include "libmesh_common.h"
#include "point.h"
#include "fem_context.h"
#include "fe_base.h"
#include "quadrature.h"
#include "elem.h"

// Local includes
#include "L-shaped.h"

#define optassert(X) {if (!(X)) libmesh_error();}

// Bring in everything from the libMesh namespace
using namespace libMesh;

void LaplaceSystem::element_qoi_derivative (DiffContext &context,
                                            const QoISet &qois)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[0]->get_JxW();

  // The basis functions for the element
  const std::vector<std::vector<Real> >          &phi = c.element_fe_var[0]->get_phi();
  
  // The element quadrature points
  const std::vector<Point > &q_point = c.element_fe_var[0]->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[0].size();
  unsigned int n_qpoints = c.element_qrule->n_points();  
  
  // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
  // we fill in the [0][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [0][0] subderivative
  DenseSubVector<Number> &Q = *c.elem_qoi_subderivatives[0][0];
      
  // Loop over the qps
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);
                  
      // If in the sub-domain over which QoI 0 is supported, add contributions
      // to the adjoint rhs
      if(fabs(x - 0.875) <= 0.125 && fabs(y - 0.125) <= 0.125)
      	{  	        	 	  
	  for (unsigned int i=0; i != n_T_dofs; i++)
	    Q(i) += JxW[qp] *phi[i][qp] ;
      	}
            
    } // end of the quadrature point qp-loop
}
