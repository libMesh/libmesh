#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem_assembly.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/boundary_info.h"

// rbOOmit includes
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_theta.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::boundary_id_type;
using libMesh::DenseSubMatrix;
using libMesh::ElemAssembly;
using libMesh::FEInterface;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBAssemblyExpansion;
using libMesh::RBConstruction;
using libMesh::RBParameters;
using libMesh::RBParametrizedFunction;
using libMesh::RBTheta;
using libMesh::RBThetaExpansion;
using libMesh::RBEIMAssembly;
using libMesh::RBEIMEvaluation;
using libMesh::RBEIMTheta;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::Elem;
using libMesh::FEBase;
using libMesh::subdomain_id_type;

// The function we're approximating with EIM
struct Gxyz : public RBParametrizedFunction
{
  unsigned int get_n_components() const
  {
    return 3;
  }

  virtual Number evaluate(const RBParameters & mu,
                          unsigned int comp,
                          const Point & p,
                          subdomain_id_type /*subdomain_id*/,
                          const std::vector<Point> & /*p_perturb*/)
  {
    Real curvature = mu.get_value("curvature");

    if(comp == 0)
    {
      return 1. + curvature*p(0);
    }
    else if(comp == 1)
    {
      return 1. + curvature*p(0);
    }
    else if(comp == 2)
    {
      return 1./(1. + curvature*p(0));
    }
    else
    {
      libmesh_error_msg("Error: Invalid comp");
    }

    return 0.;
  }
};

struct ThetaA0 : RBTheta
{
  virtual Number evaluate(const RBParameters & mu)
  {
    return mu.get_value("kappa") * mu.get_value("Bi");
  }
};

struct AssemblyA0 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext & c)
  {
    std::vector<boundary_id_type> bc_ids;
    c.get_system().get_mesh().get_boundary_info().boundary_ids (&c.get_elem(), c.side, bc_ids);
    for (std::vector<boundary_id_type>::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
      if (*b == 1 || *b == 2 || *b == 3 || *b == 4)
        {
          const unsigned int u_var = 0;

          FEBase * side_fe = nullptr;
          c.get_side_fe(u_var, side_fe);

          const std::vector<Real> & JxW_side = side_fe->get_JxW();

          const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

          // The number of local degrees of freedom in each variable
          const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

          // Now we will build the affine operator
          unsigned int n_sidepoints = c.get_side_qrule().n_points();

          for (unsigned int qp=0; qp != n_sidepoints; qp++)
            for (unsigned int i=0; i != n_u_dofs; i++)
              for (unsigned int j=0; j != n_u_dofs; j++)
                c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];

          break;
        }
  }
};

struct ThetaA1 : RBTheta
{
  virtual Number evaluate(const RBParameters & mu)
  {
    return mu.get_value("kappa") * mu.get_value("Bi") * mu.get_value("curvature");
  }
};

struct AssemblyA1 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext & c)
  {
    std::vector<boundary_id_type> bc_ids;
    c.get_system().get_mesh().get_boundary_info().boundary_ids (&c.get_elem(), c.side, bc_ids);
    for (std::vector<boundary_id_type>::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
      if (*b == 1 || *b == 3) // y == -0.2, y == 0.2
        {
          const unsigned int u_var = 0;

          FEBase * side_fe = nullptr;
          c.get_side_fe(u_var, side_fe);

          const std::vector<Real> & JxW_side = side_fe->get_JxW();

          const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

          const std::vector<Point> & xyz = side_fe->get_xyz();

          // The number of local degrees of freedom in each variable
          const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

          // Now we will build the affine operator
          unsigned int n_sidepoints = c.get_side_qrule().n_points();

          for (unsigned int qp=0; qp != n_sidepoints; qp++)
            {
              Real x_hat = xyz[qp](0);

              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW_side[qp] * x_hat * phi_side[j][qp]*phi_side[i][qp];
            }

          break;
        }
  }
};

struct ThetaA2 : RBTheta {
  virtual Number evaluate(const RBParameters & mu)
  {
    return 0.2*mu.get_value("kappa") * mu.get_value("Bi") * mu.get_value("curvature");
  }
};
struct AssemblyA2 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext & c)
  {
    std::vector<boundary_id_type> bc_ids;
    c.get_system().get_mesh().get_boundary_info().boundary_ids (&c.get_elem(), c.side, bc_ids);
    for (std::vector<boundary_id_type>::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
      if (*b == 2 || *b == 4) // x == 0.2, x == -0.2
        {
          const unsigned int u_var = 0;

          FEBase * side_fe = nullptr;
          c.get_side_fe(u_var, side_fe);

          const std::vector<Real> & JxW_side = side_fe->get_JxW();

          const std::vector<std::vector<Real>> & phi_side = side_fe->get_phi();

          // The number of local degrees of freedom in each variable
          const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

          // Now we will build the affine operator
          unsigned int n_sidepoints = c.get_side_qrule().n_points();

          if (*b==2)
            {
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                {
                  for (unsigned int i=0; i != n_u_dofs; i++)
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
                }
            }

          if (*b==4)
            {
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                {
                  for (unsigned int i=0; i != n_u_dofs; i++)
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      c.get_elem_jacobian()(i,j) -= JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
                }
            }
        }
  }
};

struct ThetaEIM : RBEIMTheta
{
  ThetaEIM(RBEIMEvaluation & rb_eim_eval_in, unsigned int index_in) :
    RBEIMTheta(rb_eim_eval_in, index_in)
  {}

  virtual Number evaluate(const RBParameters & mu)
  {
    return mu.get_value("kappa") * RBEIMTheta::evaluate(mu);
  }
};

struct AssemblyEIM : RBEIMAssembly
{
  AssemblyEIM(RBEIMEvaluation & rb_eim_eval_in,
              unsigned int basis_function_index_in) :
    RBEIMAssembly(rb_eim_eval_in,
                  basis_function_index_in)
  {}

  virtual void interior_assembly(FEMContext & c)
  {
    // PDE variable numbers
    const unsigned int u_var = 0;

    // EIM variable numbers
    const unsigned int Gx_var = 0;
    const unsigned int Gy_var = 1;
    const unsigned int Gz_var = 2;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    std::vector<Number> eim_values_Gx;
    evaluate_basis_function(c.get_elem().id(),
                            Gx_var,
                            eim_values_Gx);

    std::vector<Number> eim_values_Gy;
    evaluate_basis_function(c.get_elem().id(),
                            Gy_var,
                            eim_values_Gy);

    std::vector<Number> eim_values_Gz;
    evaluate_basis_function(c.get_elem().id(),
                            Gz_var,
                            eim_values_Gz);

    for (unsigned int qp=0; qp != c.get_element_qrule().n_points(); qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.get_elem_jacobian()(i,j) += JxW[qp] * (eim_values_Gx[qp]*dphi[i][qp](0)*dphi[j][qp](0) +
                                                   eim_values_Gy[qp]*dphi[i][qp](1)*dphi[j][qp](1) +
                                                   eim_values_Gz[qp]*dphi[i][qp](2)*dphi[j][qp](2));
  }
};


struct ThetaF0 : RBTheta
{
  virtual Number evaluate(const RBParameters &) { return 1.; }
};

struct AssemblyF0 : ElemAssembly
{

  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.get_elem_residual()(i) += JxW[qp] * (1.*phi[i][qp]);
  }
};

struct ThetaF1 : RBTheta
{
  virtual Number evaluate(const RBParameters & mu)
  {
    return mu.get_value("curvature");
  }
};

struct AssemblyF1 : ElemAssembly
{

  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

    const std::vector<Point> & xyz = elem_fe->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Real x_hat = xyz[qp](0);

        for (unsigned int i=0; i != n_u_dofs; i++)
          c.get_elem_residual()(i) += JxW[qp] * (1.*x_hat*phi[i][qp]);
      }
  }
};

struct Ex6InnerProduct : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
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
    attach_A_theta(&theta_a1);
    attach_A_theta(&theta_a2);
    attach_F_theta(&theta_f0); // Attach the rhs theta
    attach_F_theta(&theta_f1);
  }

  // The RBTheta member variables
  ThetaA0 theta_a0;
  ThetaA1 theta_a1;
  ThetaA2 theta_a2;
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
    attach_A_assembly(&assembly_a1);
    attach_A_assembly(&assembly_a2);
    attach_F_assembly(&assembly_f0); // Attach the rhs assembly
    attach_F_assembly(&assembly_f1);
  }

  // The ElemAssembly objects
  AssemblyA0 assembly_a0;
  AssemblyA1 assembly_a1;
  AssemblyA2 assembly_a2;
  AssemblyF0 assembly_f0;
  AssemblyF1 assembly_f1;
};

#endif
