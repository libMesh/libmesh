// General libMesh includes
#include "libmesh_common.h"
#include "elem.h"
#include "fe_base.h"
#include "fem_context.h"
#include "point.h"
#include "quadrature.h"

// Local includes
#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void LaplaceSystem::side_qoi_derivative (DiffContext &context,
                                         const QoISet &qois)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.side_fe_var[0]->get_JxW();

  // The basis functions for the side
  const std::vector<std::vector<RealGradient> > &dphi = c.side_fe_var[0]->get_dphi();
  
  // The side quadrature points
  const std::vector<Point > &q_point = c.side_fe_var[0]->get_xyz();

  // Get the normal to the side at each qp
  const std::vector<Point> &face_normals = c.side_fe_var[0]->get_normals();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[0].size();
  unsigned int n_qpoints = c.side_qrule->n_points();  

  // Fill the QoI RHS corresponding to this QoI. Since this is QoI 1
  // we fill in the [1][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [1][0] subderivative
  DenseSubVector<Number> &Q = *c.elem_qoi_subderivatives[1][0];

  const Real TOL = 1.e-5;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);                  
      
      // If on the sides where the boundary QoI is supported, add contributions
      // to the adjoint rhs
      if(fabs(y - 1.0) <= TOL && x > 0.0)
  	{  	  	  	  
  	  for (unsigned int i=0; i != n_T_dofs; i++)
  	    Q(i) += JxW[qp] * (dphi[i][qp] * face_normals[qp]);  	    	    	    	    	  
  	}      
     
    } // end of the quadrature point qp-loop
}
