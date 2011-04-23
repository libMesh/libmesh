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

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::DirichletDofAssembly;
using libMesh::ElemAssembly;
using libMesh::FEInterface;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBTheta;
using libMesh::Real;
using libMesh::RealGradient;

// Functors for the parameter-dependent part of the affine decomposition of the PDE
// The RHS and outputs just require a constant value of 1, so use a default RBTheta object there
struct ThetaA0 : RBTheta { virtual Number evaluate(const std::vector<Real>& )   { return 0.05;  } };
struct ThetaA1 : RBTheta { virtual Number evaluate(const std::vector<Real>& mu) { return mu[0]; } };
struct ThetaA2 : RBTheta { virtual Number evaluate(const std::vector<Real>& mu) { return mu[1]; } };

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
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
};


struct A1 : ElemAssembly
{
  // Convection in the x-direction
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();

    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
  }
};

struct A2 : ElemAssembly
{
  // Convection in the y-direction
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[u_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[u_var]->get_phi();

    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[u_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
  }
};

struct F0 : ElemAssembly
{
  // Source term, 1 throughout the domain
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
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] );
  }
};

struct OutputAssembly : ElemAssembly
{
  OutputAssembly(Real min_x_in, Real max_x_in,
                 Real min_y_in, Real max_y_in)
                :
                min_x(min_x_in),
                max_x(max_x_in),
                min_y(min_y_in),
                max_y(max_y_in)
  {}

  // Output: Average value over the region [min_x,max_x]x[min_y,max_y]
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
    unsigned int n_qpoints = c.element_qrule->n_points();
    
    Real output_area = (max_x-min_x) * (max_y-min_y);

    Point centroid = c.elem->centroid();
    if( (min_x <= centroid(0)) && (centroid(0) <= max_x) &&
        (min_y <= centroid(1)) && (centroid(1) <= max_y) )
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        for (unsigned int i=0; i != n_u_dofs; i++)
          c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / output_area;
  }
  
  // Member variables that define the output region in 2D
  Real min_x, max_x, min_y, max_y;
};

// Build up the dirichlet_dofs_set, which stores all the Dirichlet degrees of freedom
// in this problem. In this case all boundary dofs are Dirichlet.
struct Ex23DirichletDofAssembly : DirichletDofAssembly
{
  virtual void boundary_assembly(FEMContext &c)
  {
    const unsigned int u_var = 0;

    std::vector<unsigned int> side_dofs;
    FEInterface::dofs_on_side(c.elem, c.dim, c.element_fe_var[u_var]->get_fe_type(),
                              c.side, side_dofs);

    for(unsigned int ii=0; ii<side_dofs.size(); ii++)
      dirichlet_dofs_set.insert(c.dof_indices[side_dofs[ii]]);
  }
};
