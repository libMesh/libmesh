#include "L2system.h"

#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"

using namespace libMesh;

void L2System::init_data ()
{
  this->add_variable ("u", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();
}



void L2System::init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  c.get_element_fe(0)->get_JxW();
  c.get_element_fe(0)->get_phi();
  c.get_element_fe(0)->get_xyz();

  FEMSystem::init_context(context);
}


bool L2System::element_time_derivative (bool request_jacobian,
                                        DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.get_element_fe(0)->get_JxW();

  const std::vector<std::vector<Real> > &phi = c.get_element_fe(0)->get_phi();

  const std::vector<Point> &xyz = c.get_element_fe(0)->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(0).size(); 

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> &F = c.get_elem_residual(0);

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Number u = c.interior_value(0, qp);

      Number ufunc = (*goal_func)(c, xyz[qp]);

      for (unsigned int i=0; i != n_u_dofs; i++)
        F(i) += JxW[qp] * ((u - ufunc) * phi[i][qp]);
      if (request_jacobian)
        {
          const Number JxWxD = JxW[qp] *
            context.get_elem_solution_derivative();

          for (unsigned int i=0; i != n_u_dofs; i++)
            for (unsigned int j=0; j != n_u_dofs; ++j)
              K(i,j) += JxWxD * (phi[i][qp] * phi[j][qp]);
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}
