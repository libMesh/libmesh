#ifndef ASSEMBLY_H
#define ASSEMBLY_H

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"

// rbOOmit includes
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::ElemAssembly;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBAssemblyExpansion;
using libMesh::RBEIMAssembly;
using libMesh::RBEIMConstruction;
using libMesh::RBParametrizedFunction;
using libMesh::RBParameters;
using libMesh::RBTheta;
using libMesh::RBThetaExpansion;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::Elem;
using libMesh::FEBase;

struct ShiftedGaussian : public RBParametrizedFunction
{
  unsigned int get_n_components() const
  {
    return 1;
  }

  Number evaluate(const RBParameters & mu,
                  const Point & p,
                  subdomain_id_type /*subdomain_id*/)
  {
    Real center_x = mu.get_value("center_x");
    Real center_y = mu.get_value("center_y");
    return exp(-2.*(pow(center_x-p(0),2.) + pow(center_y-p(1),2.)));
  }
};

// Expansion of the PDE operator
struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters &) { return 0.05;  } };

struct A0 : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
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


struct EIM_IP_assembly : ElemAssembly
{
  // Use the L2 norm to find the best fit
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.get_elem_jacobian()(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
  }
};

struct EIM_F : RBEIMAssembly
{
  EIM_F(RBEIMConstruction & rb_eim_con_in,
        unsigned int basis_function_index_in) :
    RBEIMAssembly(rb_eim_con_in,
                  basis_function_index_in)
  {}

  virtual void interior_assembly(FEMContext & c)
  {
    // PDE variable number
    const unsigned int u_var = 0;

    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    // EIM variable number
    const unsigned int eim_var = 0;

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    std::vector<Number> eim_values;
    evaluate_basis_function(eim_var,
                            c.get_elem(),
                            c.get_element_qrule().get_points(),
                            eim_values);

    for (unsigned int qp=0; qp != c.get_element_qrule().n_points(); qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.get_elem_residual()(i) += JxW[qp] * (eim_values[qp]*phi[i][qp]);
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
