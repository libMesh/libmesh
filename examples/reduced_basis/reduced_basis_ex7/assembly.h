#ifndef __assembly_h__
#define __assembly_h__

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem_assembly.h"
#include "libmesh/quadrature_gauss.h"

// rbOOmit includes
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#define damping_epsilon 0.001
#define R_rad 12.0

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::ElemAssembly;
using libMesh::FEInterface;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBAssemblyExpansion;
using libMesh::RBParameters;
using libMesh::RBTheta;
using libMesh::RBThetaExpansion;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::MeshBase;
using libMesh::libmesh_conj;

// Functors for the parameter-dependent part of the affine decomposition of the PDE
struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(1., mu.get_value("frequency")*damping_epsilon); return val; } };
struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(-mu.get_value("frequency")*mu.get_value("frequency"), 0.); return val; } };
struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0., mu.get_value("frequency")); return val; } };
struct ThetaA3 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0.5/R_rad, mu.get_value("frequency")); return val; } };

struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0., 2.*mu.get_value("frequency")); return val; } };

struct ThetaOutput0 : RBTheta { virtual Number evaluate(const RBParameters& ) { Number val(1., 0.); return val; } };

struct AcousticsInnerProduct : ElemAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int p_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[p_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[p_var]->get_phi();

    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[p_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_p_dofs; i++)
        for (unsigned int j=0; j != n_p_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                             dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)) +
                                             (phi[j][qp]*libmesh_conj(phi[i][qp])) );
  }
};

struct A0 : ElemAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int p_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[p_var]->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi =
      c.element_fe_var[p_var]->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_p_dofs; i++)
        for (unsigned int j=0; j != n_p_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                             dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)));
  }
};


struct A1 : ElemAssembly
{
  virtual void interior_assembly(FEMContext &c)
  {
    const unsigned int p_var = 0;

    const std::vector<Real> &JxW =
      c.element_fe_var[p_var]->get_JxW();

    const std::vector<std::vector<Real> >& phi =
      c.element_fe_var[p_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

    // Now we will build the affine operator
    unsigned int n_qpoints = (c.get_element_qrule())->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_p_dofs; i++)
        for (unsigned int j=0; j != n_p_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * (phi[j][qp]*libmesh_conj(phi[i][qp]));
  }
};

struct A2 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext &c)
  {
    if( mesh->boundary_info->has_boundary_id (c.elem, c.side, 1) ) // Forcing on the horn "inlet"
    {
      const unsigned int p_var = 0;

      const std::vector<Real> &JxW_face =
        c.side_fe_var[p_var]->get_JxW();

      const std::vector<std::vector<Real> >& phi_face =
        c.side_fe_var[p_var]->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

      // Now we will build the affine operator
      unsigned int n_sidepoints = (c.side_qrule)->n_points();

      for (unsigned int qp=0; qp != n_sidepoints; qp++)
        for (unsigned int i=0; i != n_p_dofs; i++)
          for (unsigned int j=0; j != n_p_dofs; j++)
            c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
    }
  }

  MeshBase* mesh;
};

struct A3 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext &c)
  {
    if( mesh->boundary_info->has_boundary_id (c.elem, c.side, 2) ) // Radiation condition on the "bubble"
    {
      const unsigned int p_var = 0;

      const std::vector<Real> &JxW_face =
        c.side_fe_var[p_var]->get_JxW();

      const std::vector<std::vector<Real> >& phi_face =
        c.side_fe_var[p_var]->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

      // Now we will build the affine operator
      unsigned int n_sidepoints = (c.side_qrule)->n_points();

      for (unsigned int qp=0; qp != n_sidepoints; qp++)
        for (unsigned int i=0; i != n_p_dofs; i++)
          for (unsigned int j=0; j != n_p_dofs; j++)
            c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
    }
  }

  MeshBase* mesh;
};

struct F0 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext &c)
  {
    if( mesh->boundary_info->has_boundary_id (c.elem, c.side, 1) ) // Output is calculated on the horn "inlet"
    {
      const unsigned int p_var = 0;

      const std::vector<Real> &JxW_face =
        c.side_fe_var[p_var]->get_JxW();

      const std::vector<std::vector<Real> >& phi_face =
        c.side_fe_var[p_var]->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

      // Now we will build the affine operator
      unsigned int n_sidepoints = (c.side_qrule)->n_points();

      for (unsigned int qp=0; qp != n_sidepoints; qp++)
        for (unsigned int i=0; i != n_p_dofs; i++)
            c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
    }
  }

  MeshBase* mesh;
};

struct Output0 : ElemAssembly
{
  virtual void boundary_assembly(FEMContext &c)
  {
    if( mesh->boundary_info->has_boundary_id (c.elem, c.side, 1) ) // Forcing on the horn "inlet"
    {
      const unsigned int p_var = 0;

      const std::vector<Real> &JxW_face =
        c.side_fe_var[p_var]->get_JxW();

      const std::vector<std::vector<Real> >& phi_face =
        c.side_fe_var[p_var]->get_phi();

      // The number of local degrees of freedom in each variable
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

      // Now we will build the affine operator
      unsigned int n_sidepoints = (c.side_qrule)->n_points();

      for (unsigned int qp=0; qp != n_sidepoints; qp++)
        for (unsigned int i=0; i != n_p_dofs; i++)
            c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
    }
  }

  MeshBase* mesh;
};

struct AcousticsRBThetaExpansion : RBThetaExpansion
{

  /**
   * Constructor.
   */
  AcousticsRBThetaExpansion()
  {
    // set up the RBThetaExpansion object
    attach_A_theta(&theta_a_0);   // Attach the lhs theta
    attach_A_theta(&theta_a_1);
    attach_A_theta(&theta_a_2);
    attach_A_theta(&theta_a_3);

    attach_F_theta(&theta_f_0);    // Attach the rhs theta

    attach_output_theta(&theta_output_0);
  }

  // The RBTheta member variables
  ThetaA0 theta_a_0;
  ThetaA1 theta_a_1;
  ThetaA2 theta_a_2;
  ThetaA3 theta_a_3;
  ThetaF0 theta_f_0;
  ThetaOutput0 theta_output_0;
};

// Define an RBAssemblyExpansion class for this PDE
struct AcousticsRBAssemblyExpansion : RBAssemblyExpansion
{

  /**
   * Constructor.
   */
  AcousticsRBAssemblyExpansion(MeshBase& mesh_in)
  {
    A2_assembly.mesh = &mesh_in;
    A3_assembly.mesh = &mesh_in;
    F0_assembly.mesh = &mesh_in;
    Output0_assembly.mesh = &mesh_in;

    // And set up the RBAssemblyExpansion object
    attach_A_assembly(&A0_assembly); // Attach the lhs assembly
    attach_A_assembly(&A1_assembly);
    attach_A_assembly(&A2_assembly);
    attach_A_assembly(&A3_assembly);

    attach_F_assembly(&F0_assembly); // Attach the rhs assembly

    attach_output_assembly(&Output0_assembly);
  }

  // The ElemAssembly objects
  A0 A0_assembly;
  A1 A1_assembly;
  A2 A2_assembly;
  A3 A3_assembly;
  F0 F0_assembly;
  Output0 Output0_assembly;
  AcousticsInnerProduct acoustics_inner_product;
};

#endif

#endif
