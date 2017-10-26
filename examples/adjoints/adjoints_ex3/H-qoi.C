#include "H-qoi.h"

// Here we define the functions to compute the QoI (side_qoi)
// and supply the right hand side for the associated adjoint problem (side_qoi_derivative)

using namespace libMesh;

void CoupledSystemQoI::init_qoi(std::vector<Number> & sys_qoi)
{
  //Only 1 qoi to worry about
  sys_qoi.resize(1);
}

// This function supplies the right hand side for the adjoint problem
// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void CoupledSystemQoI::side_qoi_derivative (DiffContext & context,
                                            const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();

  // Side basis functions
  const std::vector<std::vector<Real>> & phi = side_fe->get_phi();
  const std::vector<Point > & q_point = side_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(1).size();

  DenseSubVector<Number> & Qu = c.get_qoi_derivatives(0, 0);
  DenseSubVector<Number> & QC = c.get_qoi_derivatives(0, 3);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_side_qrule().n_points();

  Number u = 0.;
  Number C = 0.;

  // // If side is on outlet
  if (c.has_side_boundary_id(2))
    {
      // Loop over all the qps on this side
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          Real x = q_point[qp](0);

          // If side is on left outlet
          if (x < 0.)
            {
              // Get u at the qp
              u = c.side_value(0, qp);
              C = c.side_value(3, qp);

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

// This function computes the actual QoI
void CoupledSystemQoI::side_qoi(DiffContext & context,
                                const QoISet & /* qois */)
{

  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();
  const std::vector<Point> & q_point = side_fe->get_xyz();

  // Loop over qp's, compute the function at each qp and add
  // to get the QoI
  unsigned int n_qpoints = c.get_side_qrule().n_points();

  Number dQoI_0 = 0.;
  Number u = 0.;
  Number C = 0.;

  // If side is on the left outlet
  if (c.has_side_boundary_id(2))
    {
      //Loop over all the qps on this side
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          Real x = q_point[qp](0);

          if (x < 0.)
            {
              // Get u and C at the qp
              u = c.side_value(0, qp);
              C = c.side_value(3, qp);

              dQoI_0 += JxW[qp] * -u * C;
            } // end if
        } // end quadrature loop
    } // end if on bdry

  c.get_qois()[0] += dQoI_0;
}
