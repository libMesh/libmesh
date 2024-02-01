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
#include "libmesh/utility.h"
#include "libmesh/print_trace.h"
#include "libmesh/int_range.h"
#include "libmesh/bounding_box.h"

// rbOOmit includes
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"

// C++ includes
#include <cmath>

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
using libMesh::subdomain_id_type;
using libMesh::dof_id_type;
using libMesh::Utility::pow;
using libMesh::make_range;
using libMesh::BoundingBox;

struct ShiftedGaussian : public RBParametrizedFunction
{
  unsigned int get_n_components() const override
  {
    return 1;
  }

  virtual std::vector<Number>
  evaluate(const RBParameters & mu,
           const Point & p,
           dof_id_type /*elem_id*/,
           unsigned int /*qp*/,
           subdomain_id_type /*subdomain_id*/,
           const std::vector<Point> & /*p_perturb*/,
           const std::vector<Real> & /*phi_i_qp*/) override
  {
    // // Old way, there is only 1 entry in the return vector
    // Real center_x = mu.get_value("center_x");
    // Real center_y = mu.get_value("center_y");
    // return std::vector<Number> { std::exp(-2. * (pow<2>(center_x - p(0)) + pow<2>(center_y - p(1)))) };

    // New way, there are get_n_components() * mu.n_samples() entries in
    // the return vector.  In debug mode, we verify that the same
    // number of samples are provided for all parameters when
    // RBParameters::n_samples() is called.
    std::vector<Number> ret(this->get_n_components() * mu.n_samples());
    for (auto i : make_range(mu.n_samples()))
      {
        Real center_x = mu.get_sample_value("center_x", i);
        Real center_y = mu.get_sample_value("center_y", i);
        ret[i] = std::exp(-2. * (pow<2>(center_x - p(0)) + pow<2>(center_y - p(1))));
      }
    return ret;
  }
};

// A simple Theta(mu) function which just returns a constant value
// (hence does not depend on mu). The constant must be specified when
// calling the constructor.
struct ThetaConstant : RBTheta
{
  /**
   * Constructor
   */
  ThetaConstant(Number val) : _val(val) {}

  /**
   * Evaluate theta for a single scalar-valued RBParameters object.
   * In this case, Theta(mu) does not depend on mu explicitly, except
   * to throw an error in case a multi-sample RBParameters object has
   * been provided, since the evaluate() interface does not support
   * multi-sample RBParameters objects.
   */
  virtual Number evaluate(const RBParameters & mu) override
  {
    libmesh_error_msg_if(mu.n_samples() > 1,
                         "You should only call the evaluate_vec() API when using multi-sample RBParameters objects.");

    return _val;
  }

  /**
   * Evaluate theta for multiple mu values, each of which may have
   * multiple "samples".  In this case, Theta(mu) does not depend on mu
   * explicitly, except to determine the number of "samples" (aka
   * n_samples()) which mu has, so that the output vector is sized
   * appropriately.
   */
  virtual std::vector<Number> evaluate_vec(const std::vector<RBParameters> & mus) override
  {
    // Compute the number of values to be returned in the vector. For
    // single-sample RBParameters objects, there would be mus.size()
    // values returned in the vector. For multi-sample RBParameters
    // objects, there are:
    // sum_i mus[i].n_samples()
    // total Thetas, i.e. one Theta per sample.
    unsigned int count = 0;
    for (const auto & mu : mus)
      count += mu.n_samples();

    return std::vector<Number>(count, _val);
  }

private:
  Number _val;
};

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
    evaluate_basis_function(c.get_elem().id(),
                            eim_var,
                            eim_values);

    for (unsigned int qp=0; qp != c.get_element_qrule().n_points(); qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.get_elem_residual()(i) += JxW[qp] * (eim_values[qp]*phi[i][qp]);
  }
};

/**
 * Output assembly object which computes the average value of the
 * solution variable inside a user-provided BoundingBox.
 * OutputAssembly is also used in reduced_basis_ex1.
 */
struct OutputAssembly : ElemAssembly
{
  OutputAssembly(const BoundingBox & bbox_in) :
    bbox(bbox_in)
  {}

  // Output: Average value over the region [min_x,max_x]x[min_y,max_y]
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * fe = nullptr;
    c.get_element_fe(u_var, fe);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    // TODO: BoundingBox should be able to compute and return its area/volume
    Real output_area = (bbox.max()(0) - bbox.min()(0)) * (bbox.max()(1) - bbox.min()(1));

    if (bbox.contains_point(c.get_elem().vertex_average()))
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        for (unsigned int i=0; i != n_u_dofs; i++)
          c.get_elem_residual()(i) += JxW[qp] * phi[i][qp] / output_area;
  }

  // Member variables that define the output region in 2D
  BoundingBox bbox;
};

// Define an RBThetaExpansion class for this PDE
struct EimTestRBThetaExpansion : RBThetaExpansion
{
  /**
   * Constructor.
   */
  EimTestRBThetaExpansion() :
    theta_a_0(0.05),
    output_theta(1.0)
  {
    attach_A_theta(&theta_a_0);

    // Note: there are no Thetas associated with the RHS since we use
    // an EIM approximation for the forcing term.

    // Attach an RBTheta object for each output. Here we just use the
    // ThetaConstant class again, but this time with a value of 1.0.
    attach_output_theta(&output_theta);
    attach_output_theta(&output_theta);
  }

  // The RBTheta member variables
  ThetaConstant theta_a_0;
  ThetaConstant output_theta;
};

// Define an RBAssemblyExpansion class for this PDE
struct EimTestRBAssemblyExpansion : RBAssemblyExpansion
{
  /**
   * Constructor.
   */
  EimTestRBAssemblyExpansion():
    L0(BoundingBox(/*min=*/Point(-0.2, -0.2), /*max=*/Point(0.0, 0.0))),
    L1(BoundingBox(/*min=*/Point(0.0, 0.0), /*max=*/Point(0.2, 0.2)))
  {
    attach_A_assembly(&A0_assembly);

    // Attach output assembly objects
    attach_output_assembly(&L0);
    attach_output_assembly(&L1);
  }

  // A0 assembly object
  A0 A0_assembly;

  // Assembly objects associated with the output functionals
  OutputAssembly L0;
  OutputAssembly L1;
};

#endif
