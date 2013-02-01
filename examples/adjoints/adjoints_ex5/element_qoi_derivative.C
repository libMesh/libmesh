// General libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"

// Local includes
#include "heatsystem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void HeatSystem::element_qoi_derivative (DiffContext &context,
                                            const QoISet & /* qois */)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  
  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[0]->get_JxW();

  // The basis functions for the element
  const std::vector<std::vector<Real> >          &phi = c.element_fe_var[0]->get_phi();
    
  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[0].size();
  unsigned int n_qpoints = (c.get_element_qrule())->n_points();  
  
  // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
  // we fill in the [0][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [0][0] subderivative
  DenseSubVector<Number> &Q = *c.elem_qoi_subderivatives[0][0];

  // A reference to the system context is built with
  const System & sys = c.get_system();
  
  // Get a pointer to the adjoint solution vector
  NumericVector<Number> &adjoint_solution = const_cast<System &>(sys).get_adjoint_solution(0);

  // Get the previous adjoint solution values at all the qps
  
  std::vector<Number> old_adjoint (n_qpoints, 0);
  
  c.interior_values<Real>(0, adjoint_solution, old_adjoint);

  // Our QoI depends solely on the final time, so there are no QoI contributions.
  // However, there is a contribution from the adjoint solution timestep, for the 
  // time part of the residual of the adjoint problem
  // Loop over the qps
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {                        
      for (unsigned int i=0; i != n_T_dofs; i++)
	{
	  Q(i) += -JxW[qp] * old_adjoint[qp] * phi[i][qp] ;
      	}
            
    } // end of the quadrature point qp-loop
}
