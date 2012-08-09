// local includes
#include "assembly.h"
#include "rb_classes.h"

// libMesh includes
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_vector.h"
#include "dense_subvector.h"
#include "fe.h"
#include "fe_interface.h"
#include "fe_base.h"
#include "elem_assembly.h"
#include "quadrature_gauss.h"
#include "boundary_info.h"

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
  const Real lambda_1 = E * nu / ( (1. + nu) * (1. - 2.*nu) );
  const Real lambda_2 = 0.5 * E / (1. + nu);

  return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
       + lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
}

void AssemblyA0::interior_assembly(FEMContext &c)
{
  const unsigned int n_components = rb_sys.n_vars();
  
  // make sure we have three components
  libmesh_assert(n_components == 3);
  
  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();
  
  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
  n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
  n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
  for (unsigned int C_i = 0; C_i < n_components; C_i++)
  {
    unsigned int C_j = 0;
    for (unsigned int C_k = 0; C_k < n_components; C_k++)
    {
      for (unsigned int C_l = 1; C_l < n_components; C_l++)
      {
        
        Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
        for (unsigned int qp=0; qp<n_qpoints; qp++)
        {
          for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
          {
            for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
            {
              (*c.elem_subjacobians[C_i][C_k])(i,j) += 
                JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
          }
        }
        
      }
    }
  }
  
  for (unsigned int C_i = 0; C_i < n_components; C_i++)
  {
    for (unsigned int C_j = 1; C_j < n_components; C_j++)
    {
      for (unsigned int C_k = 0; C_k < n_components; C_k++)
      {
        unsigned int C_l = 0;
          
        Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
        for (unsigned int qp=0; qp<n_qpoints; qp++)
        {
          for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
          {
            for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
            {
              (*c.elem_subjacobians[C_i][C_k])(i,j) += 
                JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
          }
        }
          
      }
    }
  }
  
}

void AssemblyA1::interior_assembly(FEMContext &c)
{
  const unsigned int n_components = rb_sys.n_vars();
  
  // make sure we have three components
  libmesh_assert(n_components == 3);
  
  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();
  
  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
  n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
  n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
  for (unsigned int C_i = 0; C_i < n_components; C_i++)
  {
    for (unsigned int C_j = 1; C_j < n_components; C_j++)
    {
      for (unsigned int C_k = 0; C_k < n_components; C_k++)
      {
        for (unsigned int C_l = 1; C_l < n_components; C_l++)
        {
          
          Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
          for (unsigned int qp=0; qp<n_qpoints; qp++)
          {
            for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
            {
              for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
              {
                (*c.elem_subjacobians[C_i][C_k])(i,j) += 
                  JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
            }
          }
          
        }
      }
    }
  }
}

void AssemblyA2::interior_assembly(FEMContext &c)
{
  const unsigned int n_components = rb_sys.n_vars();
  
  // make sure we have three components
  libmesh_assert(n_components == 3);
  
  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();
  
  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  std::vector<unsigned int> n_var_dofs(n_components);
  n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
  n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
  n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
  for (unsigned int C_i = 0; C_i < n_components; C_i++)
  {
    unsigned int C_j = 0;
    
    for (unsigned int C_k = 0; C_k < n_components; C_k++)
    {
      unsigned int C_l = 0;

      Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
      for (unsigned int qp=0; qp<n_qpoints; qp++)
      {
        for (unsigned int i=0; i<n_var_dofs[C_i]; i++)
        {
          for (unsigned int j=0; j<n_var_dofs[C_k]; j++)
          {
            (*c.elem_subjacobians[C_i][C_k])(i,j) += 
              JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
          }
        }
      }

    }
  }
}

void AssemblyF0::boundary_assembly(FEMContext &c)
{
  if(rb_sys.get_mesh().boundary_info->boundary_id(c.elem, c.side) == BOUNDARY_ID_MAX_X)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW_side =
      c.side_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi_side =
      c.side_fe_var[u_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.side_qrule->n_points();
    DenseSubVector<Number>& Fu = *c.elem_subresiduals[u_var];

    for (unsigned int qp=0; qp < n_qpoints; qp++)
      for (unsigned int i=0; i < n_u_dofs; i++)
      {
        Fu(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
      }
  }
}

void AssemblyF1::boundary_assembly(FEMContext &c)
{
  if(rb_sys.get_mesh().boundary_info->boundary_id(c.elem, c.side) == BOUNDARY_ID_MAX_X)
  {
    const unsigned int u_var = 0;
    const unsigned int v_var = 1;

    const std::vector<Real> &JxW_side =
      c.side_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi_side =
      c.side_fe_var[u_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.side_qrule->n_points();
    DenseSubVector<Number>& Fv = *c.elem_subresiduals[v_var];

    for (unsigned int qp=0; qp < n_qpoints; qp++)
      for (unsigned int i=0; i < n_v_dofs; i++)
      {
        Fv(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
      }
  }
}

void AssemblyF2::boundary_assembly(FEMContext &c)
{
  if(rb_sys.get_mesh().boundary_info->boundary_id(c.elem, c.side) == BOUNDARY_ID_MAX_X)
  {
    const unsigned int u_var = 0;
    const unsigned int w_var = 2;

    const std::vector<Real> &JxW_side =
      c.side_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi_side =
      c.side_fe_var[u_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_w_dofs = c.dof_indices_var[w_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.side_qrule->n_points();
    DenseSubVector<Number>& Fw = *c.elem_subresiduals[w_var];

    for (unsigned int qp=0; qp < n_qpoints; qp++)
      for (unsigned int i=0; i < n_w_dofs; i++)
      {
        Fw(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
      }
  }
}

void InnerProductAssembly::interior_assembly(FEMContext &c)
{
  const unsigned int u_var = rb_sys.u_var;
  const unsigned int v_var = rb_sys.v_var;
  const unsigned int w_var = rb_sys.w_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();
  const unsigned int n_w_dofs = c.dof_indices_var[w_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = (c.get_element_qrule())->n_points();
      
  DenseSubMatrix<Number>& Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number>& Kvv = *c.elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number>& Kww = *c.elem_subjacobians[w_var][w_var];
  
  for (unsigned int qp=0; qp<n_qpoints; qp++)
  {
      for (unsigned int i=0; i<n_u_dofs; i++)
        for (unsigned int j=0; j<n_u_dofs; j++)
        {
          Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        }
      
      for (unsigned int i=0; i<n_v_dofs; i++)
        for (unsigned int j=0; j<n_v_dofs; j++)
        {
          Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        }

      for (unsigned int i=0; i<n_w_dofs; i++)
        for (unsigned int j=0; j<n_w_dofs; j++)
        {
          Kww(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        }
  }
}