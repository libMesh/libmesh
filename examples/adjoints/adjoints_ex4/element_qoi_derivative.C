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

// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void LaplaceSystem::element_qoi_derivative (DiffContext & context,
                                            const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The basis functions for the element
  const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

  // The element quadrature points
  const std::vector<Point > & q_point = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
  // we fill in the [0][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [0][0] subderivative
  DenseSubVector<Number> & Q = c.get_qoi_derivatives(0, 0);

  // Loop over the qps
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      // If in the sub-domain over which QoI 0 is supported, add contributions
      // to the adjoint rhs
      if (std::abs(x - 0.875) <= 0.125 && std::abs(y - 0.125) <= 0.125)
        for (unsigned int i=0; i != n_T_dofs; i++)
          Q(i) += JxW[qp] *phi[i][qp];
    } // end of the quadrature point qp-loop
}
