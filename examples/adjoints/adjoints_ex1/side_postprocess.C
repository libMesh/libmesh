// General libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"

// Local includes
#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define the postprocess function to compute QoI 1, the integral of the the normal
// derivative of the solution over part of the boundary
void LaplaceSystem::side_postprocess(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * side_fe = nullptr;
  c.get_side_fe(0, side_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();

  const std::vector<Point > & q_point = side_fe->get_xyz();

  const std::vector<Point> & face_normals = side_fe->get_normals();

  unsigned int n_qpoints = c.get_side_qrule().n_points();

  Number dQoI_1 = 0.;

  // Loop over qp's, compute the function at each qp and add
  // to get the QoI
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      const Real TOL = 1.e-5;

      // If on the bottom horizontal bdry (y = -1)
      if (std::abs(y - 1.0) <= TOL && x > 0.0)
        {
          // Get the value of the gradient at this point
          const Gradient grad_T = c.side_gradient(0, qp);

          // Add the contribution of this qp to the integral QoI
          dQoI_1 += JxW[qp] * (grad_T * face_normals[qp]);
        }

    } // end of the quadrature point qp-loop

  computed_QoI[1] = computed_QoI[1] + dQoI_1;
}
