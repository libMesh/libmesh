#ifndef __assembly_h__
#define __assembly_h__

#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "fe_base.h"
#include "elem_assembly.h"
#include "quadrature_gauss.h"

// rbOOmit includes
#include "rb_theta.h"
#include "rb_assembly_expansion.h"
#include "rb_parametrized_function.h"
#include "rb_eim_construction.h"

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

// The first component of the function we're approximating with EIM
struct G_0 : public RBParametrizedFunction
{
  virtual Number evaluate(const RBParameters& mu,
                          const Point& p)
  {
    Real curvature = mu.get_value("curvature");
    return 1. + curvature*p(0);
  }
};

// The second component of the function we're approximating with EIM
struct G_1 : public RBParametrizedFunction
{
  virtual Number evaluate(const RBParameters& mu,
                          const Point& p)
  {
    Real curvature = mu.get_value("curvature");
    return 1./(1. + curvature*p(0));
  }
};

// Functors for the parameter-dependent part of the affine decomposition of the PDE
struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters&   ) { return 1.; } };
struct ThetaF1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("curvature"); } };

struct Ex6InnerProduct : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
};

struct Ex6EIMInnerProduct : ElemAssembly
{

  // Use the L2 norm to find the best fit
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int g0_var = 0;
    const unsigned int g1_var = 1;

    const std::vector<Real> &JxW =
      c.element_fe_var[g0_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[g0_var]->get_phi();

    const unsigned int n_u_dofs = c.dof_indices_var[g0_var].size();

    unsigned int n_qpoints = (c.get_element_qrule())->n_points();
    
    DenseSubMatrix<Number>& K00 = *c.elem_subjacobians[g0_var][g0_var];
    DenseSubMatrix<Number>& K11 = *c.elem_subjacobians[g1_var][g1_var];

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
        {
          K00(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
          K11(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
        }
  }
};

struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& ) { return 1.; } };
struct A0 : ElemAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp](1)*dphi[i][qp](1);
  }
};

struct EIM_A : RBEIMAssembly
{

  EIM_A(RBEIMConstruction& rb_eim_con_in,
        unsigned int basis_function_index_in)
  : RBEIMAssembly(rb_eim_con_in,
                  basis_function_index_in)
  {}

  virtual void interior_assembly(FEMContext &c)
  {
    // PDE variable numbers
    const unsigned int u_var = 0;
    
    // EIM variable numbers
    const unsigned int g0_var = 0;
    const unsigned int g1_var = 1;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    const std::vector<Point>& qpoints =
      c.element_fe_var[u_var]->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    std::vector<Number> eim_values_g0;
    evaluate_basis_function(g0_var,
                            *c.elem,
                            qpoints,
                            eim_values_g0);

    std::vector<Number> eim_values_g1;
    evaluate_basis_function(g1_var,
                            *c.elem,
                            qpoints,
                            eim_values_g1);

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
        {
          c.elem_jacobian(i,j) += JxW[qp] * ( eim_values_g0[qp]*dphi[i][qp](0)*dphi[j][qp](0) + 
                                              eim_values_g1[qp]*dphi[i][qp](2)*dphi[j][qp](2) );
        }
    }
  }

};

struct F0 : ElemAssembly
{

  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] );
  }
};

struct F1 : ElemAssembly
{

  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();
    
    const std::vector<Point>& xyz =
      c.element_fe_var[u_var]->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*xyz[qp](0)*phi[i][qp] );
  }
};

// Define an RBThetaExpansion class for this PDE
// The A terms depend on EIM, so we deal with them later
struct Ex6ThetaExpansion : RBThetaExpansion
{

  /**
   * Constructor.
   */
  Ex6ThetaExpansion()
  {
    attach_A_theta(&theta_a0);
    attach_F_theta(&theta_f0); // Attach the rhs theta
    attach_F_theta(&theta_f1);
  }

  // The RBTheta member variables
  ThetaA0 theta_a0;
  ThetaF0 theta_f0;
  ThetaF1 theta_f1;
};

// Define an RBAssemblyExpansion class for this PDE
// The A terms depend on EIM, so we deal with them later
struct Ex6AssemblyExpansion : RBAssemblyExpansion
{

  /**
   * Constructor.
   */
  Ex6AssemblyExpansion()
  {
    attach_A_assembly(&assembly_a0);
    attach_F_assembly(&assembly_f0); // Attach the rhs assembly
    attach_F_assembly(&assembly_f1);
  }

  // The ElemAssembly objects
  A0 assembly_a0;
  F0 assembly_f0;
  F1 assembly_f1;
};

#endif


