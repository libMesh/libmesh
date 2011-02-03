#include "L-shaped.h"
#include "fem_context.h"
#include "fe_base.h"
#include "quadrature.h"
#include "elem.h"

// Some (older) compilers do not offer full stream functionality, OStringStream works around this.

#include "o_string_stream.h"

#include <iostream>
#include "libmesh_common.h"
#include "getpot.h"

// Including files for ofstream objects

#include <fstream>
#include <sstream>

#include "point_locator_tree.h"
#include "point.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function to assemble and compute the actual QoI, needed to compute parameter sensitivities

void LaplaceSystem::element_qoi (DiffContext &context,
			      const QoISet & qois)

{  
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
            
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[0]->get_JxW();

  const std::vector<Point> &xyz = c.element_fe_var[0]->get_xyz();

  unsigned int n_qpoints = c.element_qrule->n_points();
  
  Number dQoI_0 = 0.;
  
  // Loop over quadrature points  
  
  for (unsigned int qp = 0; qp != n_qpoints; qp++)    
    {
      // Get co-ordinate locations of the current quadrature point
      const Real xf = xyz[qp](0);
      const Real yf = xyz[qp](1);

      // If in the sub-domain omega, add the contribution to the integral R
      if(fabs(xf - 0.875) <= 0.125 && fabs(yf - 0.125) <= 0.125)
      	{      
      	  // Get the solution value at the quadrature point
      	  Number T = c.interior_value(0, qp);

      	  // Update the elemental increment dR for each qp
      	  dQoI_0 += JxW[qp] * T;
      	}

    }

  // Update the computed value of the global functional R, by adding the contribution from this element

  c.elem_qoi[0] = c.elem_qoi[0] + dQoI_0;    
  
}
