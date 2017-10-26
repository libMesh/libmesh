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
void LaplaceSystem::side_qoi_derivative (DiffContext & context,
                                         const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();

  // The basis functions for the side
  const std::vector<std::vector<RealGradient>> & dphi = side_fe->get_dphi();

  // The side quadrature points
  const std::vector<Point > & q_point = side_fe->get_xyz();

  // Get the normal to the side at each qp
  const std::vector<Point> & face_normals = side_fe->get_normals();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();
  unsigned int n_qpoints = c.get_side_qrule().n_points();

  // Fill the QoI RHS corresponding to this QoI. Since this is QoI 1
  // we fill in the [1][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [1][0] subderivative
  DenseSubVector<Number> & Q = c.get_qoi_derivatives(1, 0);

  const Real TOL = 1.e-5;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      // If on the sides where the boundary QoI is supported, add contributions
      // to the adjoint rhs
      if (std::abs(y - 1.0) <= TOL && x > 0.0)
        for (unsigned int i=0; i != n_T_dofs; i++)
          Q(i) += JxW[qp] * (dphi[i][qp] * face_normals[qp]);
    } // end of the quadrature point qp-loop
}
