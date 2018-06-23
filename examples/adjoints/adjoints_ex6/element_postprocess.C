// General libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"

// Local includes
#include "poisson.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define the postprocess function to compute QoI 0, the weighted flux

void PoissonSystem::element_postprocess (DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = nullptr;
  c.get_element_fe( 0, elem_fe );

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  const std::vector<Point> & xyz = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  Number dQoI_0 = 0.;

  // Loop over quadrature points

  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Get co-ordinate locations of the current quadrature point
      const Real x = xyz[qp](0);
      const Real y = xyz[qp](1);

      Real f = -alpha * ( ( (- 4 * alpha * alpha) * exp(-alpha*x) * y * (1 - y) ) + ( -8 + ( 8 * exp(-alpha*x) ) + ( 8 * ( 1 - exp(-alpha) )* x) ) );

      Real s = x * (1 - x) * (1 - y);
      RealVectorValue U ((1-y)*(1-(2*x)),-x*(1-x));

      Gradient grad_u = c.interior_gradient(0, qp);

      // Flux with weight s = R(u^h, s) = int ( f*s - alpha*(grad_u*grad_s) ) dx
      dQoI_0 += JxW[qp] * ( (f * s) - (alpha * (U * grad_u)) );
    }

  // Update the computed value of the global functional R, by adding the contribution from this element

  computed_QoI[0] = computed_QoI[0] + dQoI_0;

}
