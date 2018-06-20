// local includes
#include "assembly.h"
#include "rb_classes.h"

// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem_assembly.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/boundary_info.h"
#include "libmesh/node.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::ElemAssembly;
using libMesh::FEInterface;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBTheta;
using libMesh::Real;
using libMesh::RealGradient;

// Kronecker delta function
inline Real kronecker_delta(unsigned int i,
                            unsigned int j)
{
  return i == j ? 1. : 0.;
}

Real elasticity_tensor(unsigned int i,
                       unsigned int j,
                       unsigned int k,
                       unsigned int l)
{
  // Define the Poisson ratio and Young's modulus
  const Real nu = 0.3;
  const Real E  = 1.;

  // Define the Lame constants (lambda_1 and lambda_2) based on nu and E
  const Real lambda_1 = E * nu / ((1. + nu) * (1. - 2.*nu));
  const Real lambda_2 = 0.5 * E / (1. + nu);

  return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l)
    + lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
}

void AssemblyA0::interior_assembly(FEMContext & c)
{
  const unsigned int n_components = rb_sys.n_vars();

  // make sure we have three components
  libmesh_assert_equal_to (n_components, 3);

  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  FEBase * elem_fe = nullptr;
  c.get_element_fe(u_var, elem_fe);

  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.get_dof_indices(u_var).size();
  n_var_dofs[v_var] = c.get_dof_indices(v_var).size();
  n_var_dofs[w_var] = c.get_dof_indices(w_var).size();

  for (unsigned int C_i = 0; C_i < n_components; C_i++)
    {
      unsigned int C_j = 0;
      for (unsigned int C_k = 0; C_k < n_components; C_k++)
        for (unsigned int C_l = 1; C_l < n_components; C_l++)
          {
            Real C_ijkl = elasticity_tensor(C_i, C_j, C_k, C_l);
            for (unsigned int qp=0; qp<n_qpoints; qp++)
              for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
                for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
                  (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                    JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
          }
    }

  for (unsigned int C_i = 0; C_i < n_components; C_i++)
    for (unsigned int C_j = 1; C_j < n_components; C_j++)
      for (unsigned int C_k = 0; C_k < n_components; C_k++)
        {
          unsigned int C_l = 0;

          Real C_ijkl = elasticity_tensor(C_i, C_j, C_k, C_l);
          for (unsigned int qp=0; qp<n_qpoints; qp++)
            for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
              for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
                (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                  JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        }

}

void AssemblyA1::interior_assembly(FEMContext & c)
{
  const unsigned int n_components = rb_sys.n_vars();

  // make sure we have three components
  libmesh_assert_equal_to (n_components, 3);

  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  FEBase * elem_fe = nullptr;
  c.get_element_fe(u_var, elem_fe);

  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.get_dof_indices(u_var).size();
  n_var_dofs[v_var] = c.get_dof_indices(v_var).size();
  n_var_dofs[w_var] = c.get_dof_indices(w_var).size();

  for (unsigned int C_i = 0; C_i < n_components; C_i++)
    for (unsigned int C_j = 1; C_j < n_components; C_j++)
      for (unsigned int C_k = 0; C_k < n_components; C_k++)
        for (unsigned int C_l = 1; C_l < n_components; C_l++)
          {
            Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
            for (unsigned int qp=0; qp<n_qpoints; qp++)
              for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
                for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
                  (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                    JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
          }
}

void AssemblyA2::interior_assembly(FEMContext & c)
{
  const unsigned int n_components = rb_sys.n_vars();

  // make sure we have three components
  libmesh_assert_equal_to (n_components, 3);

  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  FEBase * elem_fe = nullptr;
  c.get_element_fe(u_var, elem_fe);

  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.get_dof_indices(u_var).size();
  n_var_dofs[v_var] = c.get_dof_indices(v_var).size();
  n_var_dofs[w_var] = c.get_dof_indices(w_var).size();

  for (unsigned int C_i = 0; C_i < n_components; C_i++)
    {
      unsigned int C_j = 0;

      for (unsigned int C_k = 0; C_k < n_components; C_k++)
        {
          unsigned int C_l = 0;

          Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
          for (unsigned int qp=0; qp<n_qpoints; qp++)
            for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
              for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
                (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                  JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        }
    }
}

void AssemblyF0::boundary_assembly(FEMContext & c)
{
  if (rb_sys.get_mesh().get_boundary_info().has_boundary_id
      (&c.get_elem(), c.side, BOUNDARY_ID_MAX_X))
    {
      const unsigned int u_var = 0;

      FEBase * side_fe = nullptr;
      c.get_side_fe(u_var, side_fe);

      const std::vector<Real> & JxW_side = side_fe->get_JxW();

      const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

      // Now we will build the affine operator
      unsigned int n_qpoints = c.get_side_qrule().n_points();
      DenseSubVector<Number> & Fu = c.get_elem_residual(u_var);

      for (unsigned int qp=0; qp < n_qpoints; qp++)
        for (unsigned int i=0; i < n_u_dofs; i++)
          Fu(i) += JxW_side[qp] * (1. * phi_side[i][qp]);
    }
}

void AssemblyF1::boundary_assembly(FEMContext & c)
{
  if (rb_sys.get_mesh().get_boundary_info().has_boundary_id
      (&c.get_elem(), c.side, BOUNDARY_ID_MAX_X))
    {
      const unsigned int u_var = 0;
      const unsigned int v_var = 1;

      FEBase * side_fe = nullptr;
      c.get_side_fe(u_var, side_fe);

      const std::vector<Real> & JxW_side = side_fe->get_JxW();

      const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_v_dofs = c.get_dof_indices(u_var).size();

      // Now we will build the affine operator
      unsigned int n_qpoints = c.get_side_qrule().n_points();
      DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);

      for (unsigned int qp=0; qp < n_qpoints; qp++)
        for (unsigned int i=0; i < n_v_dofs; i++)
          Fv(i) += JxW_side[qp] * (1. * phi_side[i][qp]);
    }
}

void AssemblyF2::boundary_assembly(FEMContext & c)
{
  if (rb_sys.get_mesh().get_boundary_info().has_boundary_id
      (&c.get_elem(), c.side, BOUNDARY_ID_MAX_X))
    {
      const unsigned int u_var = 0;
      const unsigned int w_var = 2;

      FEBase * side_fe = nullptr;
      c.get_side_fe(u_var, side_fe);

      const std::vector<Real> & JxW_side = side_fe->get_JxW();

      const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_w_dofs = c.get_dof_indices(w_var).size();

      // Now we will build the affine operator
      unsigned int n_qpoints = c.get_side_qrule().n_points();
      DenseSubVector<Number> & Fw = c.get_elem_residual(w_var);

      for (unsigned int qp=0; qp < n_qpoints; qp++)
        for (unsigned int i=0; i < n_w_dofs; i++)
          Fw(i) += JxW_side[qp] * (1. * phi_side[i][qp]);
    }
}

void AssemblyPointLoadX::get_nodal_rhs_values(std::map<numeric_index_type, Number> & values,
                                              const System & sys,
                                              const Node & node)
{
  // First clear the values map
  values.clear();

  if (sys.get_mesh().get_boundary_info().has_boundary_id
      (&node, NODE_BOUNDARY_ID))
    {
      numeric_index_type dof_index =
        node.dof_number(sys.number(), sys.variable_number("u"), 0);
      values[dof_index] = 1.;
    }
}

void AssemblyPointLoadY::get_nodal_rhs_values(std::map<numeric_index_type, Number> & values,
                                              const System & sys,
                                              const Node & node)
{
  // First clear the values map
  values.clear();

  if (sys.get_mesh().get_boundary_info().has_boundary_id
      (&node, NODE_BOUNDARY_ID))
    {
      numeric_index_type dof_index =
        node.dof_number(sys.number(), sys.variable_number("v"), 0);
      values[dof_index] = 1.;
    }
}

void AssemblyPointLoadZ::get_nodal_rhs_values(std::map<numeric_index_type, Number> & values,
                                              const System & sys,
                                              const Node & node)
{
  // First clear the values map
  values.clear();

  if (sys.get_mesh().get_boundary_info().has_boundary_id
      (&node, NODE_BOUNDARY_ID))
    {
      numeric_index_type dof_index =
        node.dof_number(sys.number(), sys.variable_number("w"), 0);
      values[dof_index] = 1.;
    }
}

void InnerProductAssembly::interior_assembly(FEMContext & c)
{
  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  FEBase * elem_fe = nullptr;
  c.get_element_fe(u_var, elem_fe);

  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>>& dphi = elem_fe->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();
  const unsigned int n_v_dofs = c.get_dof_indices(v_var).size();
  const unsigned int n_w_dofs = c.get_dof_indices(w_var).size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(u_var, u_var);
  DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);
  DenseSubMatrix<Number> & Kww = c.get_elem_jacobian(w_var, w_var);

  for (unsigned int qp=0; qp<n_qpoints; qp++)
    {
      for (unsigned int i=0; i<n_u_dofs; i++)
        for (unsigned int j=0; j<n_u_dofs; j++)
          Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      for (unsigned int i=0; i<n_v_dofs; i++)
        for (unsigned int j=0; j<n_v_dofs; j++)
          Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      for (unsigned int i=0; i<n_w_dofs; i++)
        for (unsigned int j=0; j<n_w_dofs; j++)
          Kww(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
    }
}
