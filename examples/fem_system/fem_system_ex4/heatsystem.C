#include "heatsystem.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"
#include "libmesh/elem.h"

using namespace libMesh;

void HeatSystem::init_data ()
{
  T_var = this->add_variable("T", static_cast<Order>(_fe_order),
                             Utility::string_to_enum<FEFamily>(_fe_family));

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Add dirichlet boundaries on all but the boundary element side
  const boundary_id_type all_ids[6] = {0, 1, 2, 3, 4, 5};
  std::set<boundary_id_type> nonyplus_bdys(all_ids, all_ids+(dim*2));
  const boundary_id_type yplus_id = (dim == 3) ? 3 : 2;
  nonyplus_bdys.erase(yplus_id);

  std::vector<unsigned int> T_only(1, T_var);
  ZeroFunction<Number> zero;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (nonyplus_bdys, T_only, zero, LOCAL_VARIABLE_ORDER));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // The temperature is evolving, with a first-order time derivative
  this->time_evolving(T_var, 1);
}



void HeatSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const std::set<unsigned char> & elem_dims =
    c.elem_dimensions();

  for (std::set<unsigned char>::const_iterator dim_it =
         elem_dims.begin(); dim_it != elem_dims.end(); ++dim_it)
    {
      const unsigned char dim = *dim_it;

      FEBase * fe = nullptr;

      c.get_element_fe(T_var, fe, dim);

      fe->get_JxW();  // For integration
      fe->get_dphi(); // For bilinear form
      fe->get_xyz();  // For forcing
      fe->get_phi();  // For forcing
    }

  FEMSystem::init_context(context);
}


bool HeatSystem::element_time_derivative (bool request_jacobian,
                                          DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const unsigned int mesh_dim =
    c.get_system().get_mesh().mesh_dimension();

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  const unsigned short dim = c.get_elem().dim();
  FEBase * fe = nullptr;
  c.get_element_fe(T_var, fe, dim);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<Point> & xyz = fe->get_xyz();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.n_dof_indices(T_var);

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(T_var, T_var);
  DenseSubVector<Number> & F = c.get_elem_residual(T_var);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution gradient at the Newton iterate
      Gradient grad_T = c.interior_gradient(T_var, qp);

      const Number k = _k[dim];

      const Point & p = xyz[qp];

      // solution + laplacian depend on problem dimension
      const Number u_exact = (mesh_dim == 2) ?
        std::sin(libMesh::pi*p(0)) * std::sin(libMesh::pi*p(1)) :
        std::sin(libMesh::pi*p(0)) * std::sin(libMesh::pi*p(1)) *
        std::sin(libMesh::pi*p(2));

      // Only apply forcing to interior elements
      const Number forcing = (dim == mesh_dim) ?
        -k * u_exact * (dim * libMesh::pi * libMesh::pi) : 0;

      const Number JxWxNK = JxW[qp] * -k;

      for (unsigned int i=0; i != n_T_dofs; i++)
        F(i) += JxWxNK * (grad_T * dphi[i][qp] + forcing * phi[i][qp]);
      if (request_jacobian)
        {
          const Number JxWxNKxD = JxWxNK *
            context.get_elem_solution_derivative();

          for (unsigned int i=0; i != n_T_dofs; i++)
            for (unsigned int j=0; j != n_T_dofs; ++j)
              K(i,j) += JxWxNKxD * (dphi[i][qp] * dphi[j][qp]);
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}
