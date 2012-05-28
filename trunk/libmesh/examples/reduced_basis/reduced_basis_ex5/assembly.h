#ifndef __assembly_h__
#define __assembly_h__

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

// rbOOmit includes
#include "rb_theta.h"
#include "rb_assembly_expansion.h"

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

// Functors for the parameter-dependent part of the affine decomposition of the PDE
struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return mu.get_value("y_scaling"); } };
struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& )     { return 1.; } };
struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return 1./mu.get_value("y_scaling"); } };
struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return mu.get_value("y_scaling") * mu.get_value("x_load"); } };
struct ThetaF1 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return mu.get_value("y_scaling") * mu.get_value("y_load"); } };

// Provide a simple subclass that just provides the elasticity tensor
struct ElasticityAssembly : ElemAssembly
{

  Real eval_elasticity_tensor(unsigned int i,
                              unsigned int j,
                              unsigned int k,
                              unsigned int l)
  {
    // Define the Poisson ratio
    const Real nu = 0.3;
  
    // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
    const Real lambda_1 = nu / ( (1. + nu) * (1. - 2.*nu) );
    const Real lambda_2 = 0.5 / (1 + nu);

    // Define the Kronecker delta functions that we need here
    Real delta_ij = (i == j) ? 1. : 0.;
    Real delta_il = (i == l) ? 1. : 0.;
    Real delta_ik = (i == k) ? 1. : 0.;
    Real delta_jl = (j == l) ? 1. : 0.;
    Real delta_jk = (j == k) ? 1. : 0.;
    Real delta_kl = (k == l) ? 1. : 0.;
  
    return lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
  }

};

struct A0 : ElasticityAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;
    const unsigned int v_var = 1;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
    const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.element_qrule->n_points();
        
    DenseSubMatrix<Number>& Kuu = *c.elem_subjacobians[u_var][u_var];
    DenseSubMatrix<Number>& Kuv = *c.elem_subjacobians[u_var][v_var];
    DenseSubMatrix<Number>& Kvu = *c.elem_subjacobians[v_var][u_var];
    DenseSubMatrix<Number>& Kvv = *c.elem_subjacobians[v_var][v_var];

      for (unsigned int qp=0; qp<n_qpoints; qp++)
      {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=0;

              C_j=0, C_l=0;
              Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=1;

              C_j=0, C_l=0;
              Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=0;

              C_j=0, C_l=0;
              Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=1;

              C_j=0, C_l=0;
              Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
      }
  }
};

struct A1 : ElasticityAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;
    const unsigned int v_var = 1;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
    const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.element_qrule->n_points();
        
    DenseSubMatrix<Number>& Kuu = *c.elem_subjacobians[u_var][u_var];
    DenseSubMatrix<Number>& Kuv = *c.elem_subjacobians[u_var][v_var];
    DenseSubMatrix<Number>& Kvu = *c.elem_subjacobians[v_var][u_var];
    DenseSubMatrix<Number>& Kvv = *c.elem_subjacobians[v_var][v_var];

      for (unsigned int qp=0; qp<n_qpoints; qp++)
      {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=0;
              
              C_j=1, C_l=0;
              Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

              C_j=0, C_l=1;
              Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=1;
              
              C_j=1, C_l=0;
              Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

              C_j=0, C_l=1;
              Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=0;
              
              C_j=1, C_l=0;
              Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

              C_j=0, C_l=1;
              Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=1;
              
              C_j=1, C_l=0;
              Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

              C_j=0, C_l=1;
              Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
      }
  }
};

struct A2 : ElasticityAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;
    const unsigned int v_var = 1;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
    const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.element_qrule->n_points();
        
    DenseSubMatrix<Number>& Kuu = *c.elem_subjacobians[u_var][u_var];
    DenseSubMatrix<Number>& Kuv = *c.elem_subjacobians[u_var][v_var];
    DenseSubMatrix<Number>& Kvu = *c.elem_subjacobians[v_var][u_var];
    DenseSubMatrix<Number>& Kvv = *c.elem_subjacobians[v_var][v_var];

      for (unsigned int qp=0; qp<n_qpoints; qp++)
      {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=0;

              C_j=1, C_l=1;
              Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=0, C_k=1;

              C_j=1, C_l=1;
              Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=0;

              C_j=1, C_l=1;
              Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i, C_j, C_k, C_l;
              C_i=1, C_k=1;

              C_j=1, C_l=1;
              Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
      }
  }
};

struct InnerProductAssembly : ElemAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    A0_assembly.interior_assembly(c);
    A1_assembly.interior_assembly(c);
    A2_assembly.interior_assembly(c);
  }
  
  A0 A0_assembly;
  A1 A1_assembly;
  A2 A2_assembly;
};

struct F0 : ElemAssembly
{
  // Apply a traction 
  virtual void boundary_assembly(FEMContext &c)
  {
    if(c.side == 1)
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

      for (unsigned int qp=0; qp != n_qpoints; qp++)
        for (unsigned int i=0; i != n_u_dofs; i++)
        {
          (*c.elem_subresiduals[u_var])(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
        }
    }
  }
};

struct F1 : ElemAssembly
{
  // Apply a traction 
  virtual void boundary_assembly(FEMContext &c)
  {
    if(c.side == 1)
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

      for (unsigned int qp=0; qp != n_qpoints; qp++)
        for (unsigned int i=0; i != n_v_dofs; i++)
        {
          (*c.elem_subresiduals[v_var])(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
        }
    }
  }
};

// Define an RBThetaExpansion class for this PDE
struct ElasticityThetaExpansion : RBThetaExpansion
{

  /**
   * Constructor.
   */
  ElasticityThetaExpansion()
  {
    // set up the RBThetaExpansion object
    attach_theta_q_a(&theta_a_0);
    attach_theta_q_a(&theta_a_1);
    attach_theta_q_a(&theta_a_2);
    attach_theta_q_f(&theta_f_0);
    attach_theta_q_f(&theta_f_1);
  }

  // The RBTheta member variables
  ThetaA0 theta_a_0;
  ThetaA1 theta_a_1;
  ThetaA2 theta_a_2;
  ThetaF0 theta_f_0;
  ThetaF1 theta_f_1;
};

// Define an RBAssemblyExpansion class for this PDE
struct ElasticityAssemblyExpansion : RBAssemblyExpansion
{

  /**
   * Constructor.
   */
  ElasticityAssemblyExpansion()
  {
    // And set up the RBAssemblyExpansion object
    attach_A_q_assembly(&A0_assembly);
    attach_A_q_assembly(&A1_assembly);
    attach_A_q_assembly(&A2_assembly);
    attach_F_q_assembly(&F0_assembly);
    attach_F_q_assembly(&F1_assembly);
  }

  // The ElemAssembly objects
  A0 A0_assembly;
  A1 A1_assembly;
  A2 A2_assembly;
  F0 F0_assembly;
  F1 F1_assembly;
};

#endif


