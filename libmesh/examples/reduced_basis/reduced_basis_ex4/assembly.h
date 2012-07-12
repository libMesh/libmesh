#ifndef __assembly_h__
#define __assembly_h__

// libMesh includes
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature.h"
#include "dof_map.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "fe_interface.h"
#include "elem.h"

// rbOOmit includes
#include "rb_assembly_expansion.h"
#include "rb_eim_theta.h"
#include "rb_parametrized_function.h"

struct ShiftedGaussian : public RBParametrizedFunction
{
  virtual Number evaluate(const RBParameters& mu,
                          const Point& p)
  {
    Real center_x = mu.get_value("center_x");
    Real center_y = mu.get_value("center_y");
    return exp( -2.*(pow(center_x-p(0),2.) + pow(center_y-p(1),2.)) );
  }
};

// Expansion of the PDE operator
struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& )   { return 0.05;  } };

struct A0 : ElemAssembly
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


struct EIM_IP_assembly : ElemAssembly
{

  // Use the L2 norm to find the best fit
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();

    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
  }
};

struct EIM_F : RBEIMAssembly
{

  EIM_F(RBEIMConstruction& rb_eim_con_in,
        unsigned int basis_function_index_in)
  : RBEIMAssembly(rb_eim_con_in,
                  basis_function_index_in)
  {}

  virtual void interior_assembly(FEMContext &c)
  {
    // PDE variable number
    const unsigned int u_var = 0;
    
    // EIM variable number
    const unsigned int eim_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();
    
    const std::vector<Point>& qpoints =
      c.element_fe_var[u_var]->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    std::vector<Number> eim_values;
    evaluate_basis_function(eim_var,
                            *c.elem,
                            qpoints,
                            eim_values);

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( eim_values[qp]*phi[i][qp] );
  }

};

// Define an RBThetaExpansion class for this PDE
struct EimTestRBThetaExpansion : RBThetaExpansion
{

  /**
   * Constructor.
   */
  EimTestRBThetaExpansion()
  {
    attach_A_theta(&theta_a_0);
  }

  // The RBTheta member variables
  ThetaA0 theta_a_0;
};

// Define an RBAssemblyExpansion class for this PDE
struct EimTestRBAssemblyExpansion : RBAssemblyExpansion
{

  /**
   * Constructor.
   */
  EimTestRBAssemblyExpansion()
  {
    attach_A_assembly(&A0_assembly);
  }

  // A0 assembly object
  A0 A0_assembly;

};

#endif
