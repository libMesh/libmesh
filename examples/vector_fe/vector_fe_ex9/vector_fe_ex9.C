// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Vector Finite Elements Example 8 - Hybridizable Discontinuous Galerkin Navier Stokes</h1>
// \author Alexander Lindsay
// \date 2023

// Basic utilities.
#include "libmesh/string_to_enum.h"

// The solver packages supported by libMesh.
#include "libmesh/enum_solver_package.h"

// The mesh object and mesh generation and modification utilities.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"

// Matrix and vector types.
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"

// The finite element object and the geometric element type.
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"

// Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "libmesh/enum_quadrature_type.h"

// The dof map, which handles degree of freedom indexing.
#include "libmesh/dof_map.h"

// The system of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

// The exact solution and error computation.
#include "libmesh/exact_solution.h"
#include "libmesh/enum_norm_type.h"
#include "solution_function.h"

// I/O utilities.
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"

#ifdef LIBMESH_HAVE_EIGEN_DENSE
#ifdef LIBMESH_HAVE_METAPHYSICL

#include <Eigen/Dense>

using namespace libMesh;
using namespace Eigen;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
typedef MatrixXcd EigenMatrix;
typedef VectorXcd EigenVector;
#else
typedef MatrixXd EigenMatrix;
typedef VectorXd EigenVector;
#endif

namespace libMesh
{
class NonlinearImplicitSystem;
}

class USoln
{
public:
  USoln(const Real mu_in) : mu(mu_in) {}

  Real operator()(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return sin(1. / 2 * y * pi) * cos(1. / 2 * x * pi);
  }

  Real forcing(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return (1. / 2) * pi * pi * mu * sin((1. / 2) * y * pi) * cos((1. / 2) * x * pi) -
           1. / 4 * pi * sin((1. / 4) * x * pi) * sin((3. / 2) * y * pi);
  }

private:
  const Real mu;
};

class VSoln
{
public:
  VSoln(const Real mu_in) : mu(mu_in) {}

  Real operator()(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return sin((1. / 4) * x * pi) * cos((1. / 2) * y * pi);
  }

  Real forcing(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return (5. / 16) * pi * pi * mu * sin((1. / 4) * x * pi) * cos((1. / 2) * y * pi) +
           (3. / 2) * pi * cos((1. / 4) * x * pi) * cos((3. / 2) * y * pi);
  }

private:
  const Real mu;
};

class PSoln
{
public:
  PSoln() = default;

  Real operator()(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return sin((3. / 2) * y * pi) * cos((1. / 4) * x * pi);
  }

  Real forcing(const Point & p) const
  {
    const auto x = p(0);
    const auto y = p(1);
    return -1. / 2 * pi * sin((1. / 4) * x * pi) * sin((1. / 2) * y * pi) -
           1. / 2 * pi * sin((1. / 2) * x * pi) * sin((1. / 2) * y * pi);
  }
};

// compute a solution indexable at quadrature points composed from the local degree of freedom
// solution vector and associated basis functions
template <typename SolnType, typename PhiType>
void
compute_qp_soln(std::vector<SolnType> & qp_vec,
                const unsigned int n_qps,
                const std::vector<std::vector<PhiType>> & phi,
                const std::vector<Number> & dof_values)
{
  libmesh_assert(dof_values.size() == phi.size());
  qp_vec.resize(n_qps);
  for (auto & val : qp_vec)
    val = {};
  for (const auto i : index_range(phi))
  {
    const auto & qp_phis = phi[i];
    libmesh_assert(qp_phis.size() == n_qps);
    const auto sol = dof_values[i];
    for (const auto qp : make_range(n_qps))
      qp_vec[qp] += qp_phis[qp] * sol;
  }
}

class HDGProblem : public libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian,
                   public libMesh::NonlinearImplicitSystem::ComputePostCheck
{
public:
  HDGProblem() : u_true_soln(mu), v_true_soln(mu) {}

  System * mixed_system;
  ImplicitSystem * lm_system;
  const MeshBase * mesh;
  const DofMap * mixed_dof_map;
  const DofMap * lm_dof_map;
  std::unique_ptr<FEVectorBase> vector_fe;
  std::unique_ptr<FEBase> scalar_fe;
  std::unique_ptr<QBase> qrule;
  std::unique_ptr<FEVectorBase> vector_fe_face;
  std::unique_ptr<FEBase> scalar_fe_face;
  std::unique_ptr<FEBase> lm_fe_face;
  std::unique_ptr<QBase> qface;
  // Parallel version of LM increment
  NumericVector<Number> * parallel_increment;
  // Ghosted version of LM increment
  NumericVector<Number> * ghosted_increment;
  // Ghosted version of old solution
  NumericVector<Number> * ghosted_old_solution;
  SparseMatrix<Number> * J;
  NumericVector<Number> * residual;
  boundary_id_type left_bnd;
  boundary_id_type top_bnd;
  boundary_id_type right_bnd;
  boundary_id_type bottom_bnd;

  void init()
  {
    // Attach quadrature rules for the FE objects that we will reinit within the element "volume"
    vector_fe->attach_quadrature_rule(qrule.get());
    scalar_fe->attach_quadrature_rule(qrule.get());

    // Attach quadrature rules for the FE objects that we will reinit on the element faces
    vector_fe_face->attach_quadrature_rule(qface.get());
    scalar_fe_face->attach_quadrature_rule(qface.get());
    lm_fe_face->attach_quadrature_rule(qface.get());

    // pre-request our required volumetric data
    JxW = &vector_fe->get_JxW();
    q_point = &vector_fe->get_xyz();
    vector_phi = &vector_fe->get_phi();
    scalar_phi = &scalar_fe->get_phi();
    grad_scalar_phi = &scalar_fe->get_dphi();
    div_vector_phi = &vector_fe->get_div_phi();

    // pre-request our required element face data
    vector_phi_face = &vector_fe_face->get_phi();
    scalar_phi_face = &scalar_fe_face->get_phi();
    lm_phi_face = &lm_fe_face->get_phi();
    JxW_face = &vector_fe_face->get_JxW();
    qface_point = &vector_fe_face->get_xyz();
    normals = &vector_fe_face->get_normals();

    const auto & bnd_info = mesh->get_boundary_info();
    left_bnd = bnd_info.get_id_by_name("left");
    top_bnd = bnd_info.get_id_by_name("top");
    right_bnd = bnd_info.get_id_by_name("right");
    bottom_bnd = bnd_info.get_id_by_name("bottom");
    libmesh_assert(left_bnd != BoundaryInfo::invalid_id);
    libmesh_assert(top_bnd != BoundaryInfo::invalid_id);
    libmesh_assert(right_bnd != BoundaryInfo::invalid_id);
    libmesh_assert(bottom_bnd != BoundaryInfo::invalid_id);
  }

  virtual void residual_and_jacobian(const NumericVector<Number> & /*X*/,
                                     NumericVector<Number> * R_in,
                                     SparseMatrix<Number> * J_in,
                                     NonlinearImplicitSystem & /*S*/) override
  {
    residual = R_in;
    J = J_in;
    assemble(true);
  }

  virtual void postcheck(const NumericVector<Number> & old_soln,
                         NumericVector<Number> & /*search_direction*/,
                         NumericVector<Number> & new_soln,
                         bool & /*changed_search_direction*/,
                         bool & /*changed_new_soln*/,
                         NonlinearImplicitSystem & /*S*/) override
  {
    *parallel_increment = new_soln;
    *parallel_increment -= old_soln;
    *ghosted_increment = *parallel_increment;
    *ghosted_old_solution = old_soln;
    assemble(false);
  }

private:
  void create_identity_residual(const QBase & quadrature,
                                const std::vector<Real> & JxW,
                                const std::vector<std::vector<Real>> & phi,
                                const std::vector<Number> & sol,
                                const std::size_t n_dofs,
                                const unsigned int i_offset)
  {
    for (const auto qp : make_range(quadrature.n_points()))
      for (const auto i : make_range(n_dofs))
        LMVec(i_offset + i) -= JxW[qp] * phi[i][qp] * sol[qp];
  }

  void create_identity_jacobian(const QBase & quadrature,
                                const std::vector<Real> & JxW,
                                const std::vector<std::vector<Real>> & phi,
                                const std::size_t n_dofs,
                                const unsigned int ij_offset)
  {
    for (const auto qp : make_range(quadrature.n_points()))
      for (const auto i : make_range(n_dofs))
        for (const auto j : make_range(n_dofs))
          LMMat(ij_offset + i, ij_offset + j) -= JxW[qp] * phi[i][qp] * phi[j][qp];
  }

  void compute_stress(const std::vector<Gradient> & vel_gradient,
                      const std::vector<Number> & p_sol,
                      const unsigned int vel_component,
                      std::vector<Gradient> & sigma)
  {
    sigma.resize(qrule->n_points());
    for (const auto qp : make_range(qrule->n_points()))
    {
      Gradient qp_p;
      qp_p(vel_component) = p_sol[qp];
      sigma[qp] = mu * vel_gradient[qp] - qp_p;
    }
  }

  void vector_volume_residual(const unsigned int i_offset,
                              const std::vector<Gradient> & vector_sol,
                              const std::vector<Number> & scalar_sol)
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(vector_n_dofs))
      {
        // Vector equation dependence on vector dofs
        MixedVec(i_offset + i) += (*JxW)[qp] * ((*vector_phi)[i][qp] * vector_sol[qp]);

        // Vector equation dependence on scalar dofs
        MixedVec(i_offset + i) += (*JxW)[qp] * ((*div_vector_phi)[i][qp] * scalar_sol[qp]);
      }
  }

  void vector_volume_jacobian(const unsigned int i_offset,
                              const unsigned int vector_j_offset,
                              const unsigned int scalar_j_offset)
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(vector_n_dofs))
      {
        // Vector equation dependence on vector dofs
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(i_offset + i, vector_j_offset + j) +=
              (*JxW)[qp] * ((*vector_phi)[i][qp] * (*vector_phi)[j][qp]);

        // Vector equation dependence on scalar dofs
        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(i_offset + i, scalar_j_offset + j) +=
              (*JxW)[qp] * ((*div_vector_phi)[i][qp] * (*scalar_phi)[j][qp]);
      }
  }

  void scalar_volume_residual(const unsigned int i_offset,
                              const std::vector<Gradient> & vel_gradient,
                              const unsigned int vel_component,
                              std::vector<Gradient> & sigma)
  {
    // const auto dim = mesh->mesh_dimension();
    compute_stress(vel_gradient, p_sol, vel_component, sigma);
    for (const auto qp : make_range(qrule->n_points()))
    {
      // // Prepare forcing function
      // const Real x = (*q_point)[qp](0);
      // const Real y = (*q_point)[qp](1);
      // const Real z = (*q_point)[qp](2);

      // // "f" is the forcing function for the Poisson equation, which is
      // // just the divergence of the exact solution for the vector field.
      // // This is the well-known "method of manufactured solutions".
      // Real f = 0;
      // if (dim == 2)
      //   f = MixedExactSolution().forcing(x, y);
      // else if (dim == 3)
      //   f = MixedExactSolution().forcing(x, y, z);

      for (const auto i : make_range(scalar_n_dofs))
      {
        // Scalar equation dependence on vector dofs
        MixedVec(i_offset + i) += (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * sigma[qp]);

        // // Scalar equation RHS
        // MixedVec(i_offset + i) += (*JxW)[qp] * (*scalar_phi)[i][qp] * f;
      }
    }
  }

  void scalar_volume_jacobian(const unsigned int i_offset,
                              const unsigned int vel_gradient_j_offset,
                              const unsigned int p_j_offset,
                              const unsigned int vel_component)
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        // Scalar equation dependence on vector dofs
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(i_offset + i, vel_gradient_j_offset + j) +=
              (*JxW)[qp] * mu * ((*grad_scalar_phi)[i][qp] * (*vector_phi)[j][qp]);

        // Scalar equation dependence on pressure dofs
        for (const auto j : make_range(p_n_dofs))
        {
          Gradient p_phi;
          p_phi(vel_component) = (*scalar_phi)[j][qp];
          MixedLM(i_offset + i, p_j_offset + j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * p_phi);
        }
      }
  }

  void pressure_volume_residual(const unsigned int i_offset,
                                const std::vector<Number> & u_sol,
                                const std::vector<Number> & v_sol)
  {
    for (const auto qp : make_range(qrule->n_points()))
    {
      const Gradient vel(u_sol[qp], v_sol[qp]);
      for (const auto i : make_range(p_n_dofs))
        LMVec(i_offset + i) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vel);
    }
  }

  void pressure_volume_jacobian(const unsigned int i_offset,
                                const unsigned int u_j_offset,
                                const unsigned int v_j_offset)
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(p_n_dofs))
        for (const auto j : make_range(scalar_n_dofs))
        {
          {
            const Gradient phi((*scalar_phi)[j][qp], 0);
            LMMixed(i_offset + i, u_j_offset + j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * phi);
          }
          {
            const Gradient phi(0, (*scalar_phi)[j][qp]);
            LMMixed(i_offset + i, v_j_offset + j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * phi);
          }
        }
  }

  void pressure_face_residual(const unsigned int i_offset,
                              const std::vector<Number> & lm_u_sol,
                              const std::vector<Number> & lm_v_sol)
  {
    for (const auto qp : make_range(qface->n_points()))
    {
      const Gradient vel(lm_u_sol[qp], lm_v_sol[qp]);
      const auto vdotn = vel * (*normals)[qp];
      for (const auto i : make_range(p_n_dofs))
        LMVec(i_offset + i) += vdotn * (*scalar_phi_face)[i][qp];
    }
  }

  void pressure_face_jacobian(const unsigned int i_offset,
                              const unsigned int lm_u_j_offset,
                              const unsigned int lm_v_j_offset)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(p_n_dofs))
        for (const auto j : make_range(lm_n_dofs))
        {
          {
            const Gradient phi((*lm_phi_face)[j][qp], 0);
            LMMat(i_offset + i, lm_u_j_offset + j) +=
                phi * (*normals)[qp] * (*scalar_phi_face)[i][qp];
          }
          {
            const Gradient phi(0, (*lm_phi_face)[j][qp]);
            LMMat(i_offset + i, lm_v_j_offset + j) +=
                phi * (*normals)[qp] * (*scalar_phi_face)[i][qp];
          }
        }
  }

  void pressure_dirichlet_residual(const unsigned int i_offset,
                                   const RealVectorValue & dirichlet_velocity)
  {
    for (const auto qp : make_range(qface->n_points()))
    {
      const auto vdotn = dirichlet_velocity * (*normals)[qp];
      for (const auto i : make_range(p_n_dofs))
        LMVec(i_offset + i) += vdotn * (*scalar_phi_face)[i][qp];
    }
  }

  void vector_dirichlet_residual(const unsigned int i_offset, const Real scalar_value)
  {
    for (const auto qp : make_range(qface->n_points()))
      // External boundary -> Dirichlet faces -> Vector equation RHS
      for (const auto i : make_range(vector_n_dofs))
        MixedVec(i_offset + i) -=
            (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * scalar_value;
  }

  void vector_face_residual(const unsigned int i_offset, const std::vector<Number> & lm_sol)
  {
    for (const auto qp : make_range(qface->n_points()))
      // Vector equation dependence on LM dofs
      for (const auto i : make_range(vector_n_dofs))
        MixedVec(i_offset + i) -=
            (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * lm_sol[qp];
  }

  void vector_face_jacobian(const unsigned int i_offset, const unsigned int lm_j_offset)
  {
    for (const auto qp : make_range(qface->n_points()))
      // Vector equation dependence on LM dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(lm_n_dofs))
          MixedLM(i_offset + i, lm_j_offset + j) -= (*JxW_face)[qp] *
                                                    ((*vector_phi_face)[i][qp] * (*normals)[qp]) *
                                                    (*lm_phi_face)[j][qp];
  }

  void scalar_dirichlet_residual(const unsigned int i_offset,
                                 const std::vector<Gradient> & vector_sol,
                                 const std::vector<Number> & scalar_sol,
                                 const unsigned int vel_component,
                                 const Real scalar_value)
  {
    for (const auto qp : make_range(qface->n_points()))
    {
      Gradient qp_p;
      qp_p(vel_component) = p_sol[qp];

      for (const auto i : make_range(scalar_n_dofs))
      {
        // vector
        MixedVec(i_offset + i) -=
            (*JxW_face)[qp] * mu * (*scalar_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // pressure
        MixedVec(i_offset + i) +=
            (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

        // scalar
        MixedVec(i_offset + i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau *
                                  scalar_sol[qp] * (*normals)[qp] * (*normals)[qp];

        // dirichlet
        MixedVec(i_offset + i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * scalar_value *
                                  (*normals)[qp] * (*normals)[qp];
      }
    }
  }

  void scalar_dirichlet_jacobian(const unsigned int i_offset,
                                 const unsigned int vector_j_offset,
                                 const unsigned int scalar_j_offset,
                                 const unsigned int p_j_offset,
                                 const unsigned int vel_component)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(i_offset + i, vector_j_offset + j) -=
              (*JxW_face)[qp] * mu * (*scalar_phi_face)[i][qp] *
              ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(p_n_dofs))
        {
          Gradient p_phi;
          p_phi(vel_component) = (*scalar_phi_face)[j][qp];
          // pressure
          MixedLM(i_offset + i, p_j_offset + j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
        }

        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(i_offset + i, scalar_j_offset + j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
              (*normals)[qp] * (*normals)[qp];
      }
  }

  void scalar_outlet_residual(const unsigned int i_offset,
                              const std::vector<Number> & scalar_sol,
                              const std::vector<Number> & lm_sol)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        // scalar
        MixedVec(i_offset + i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau *
                                  scalar_sol[qp] * (*normals)[qp] * (*normals)[qp];

        // lm
        MixedVec(i_offset + i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * lm_sol[qp] *
                                  (*normals)[qp] * (*normals)[qp];
      }
  }

  void scalar_outlet_jacobian(const unsigned int i_offset,
                              const unsigned int scalar_j_offset,
                              const unsigned int lm_j_offset)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(i_offset + i, scalar_j_offset + j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
              (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          MixedLM(i_offset + i, lm_j_offset + j) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] *
                                                    tau * (*lm_phi_face)[j][qp] * (*normals)[qp] *
                                                    (*normals)[qp];
      }
  }

  void scalar_face_residual(const unsigned int i_offset,
                            const std::vector<Gradient> & vector_sol,
                            const std::vector<Number> & scalar_sol,
                            const std::vector<Number> & lm_sol,
                            const unsigned int vel_component)
  {
    for (const auto qp : make_range(qface->n_points()))
    {
      Gradient qp_p;
      qp_p(vel_component) = p_sol[qp];
      for (const auto i : make_range(scalar_n_dofs))
      {
        // vector
        MixedVec(i_offset + i) -=
            (*JxW_face)[qp] * mu * (*scalar_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // pressure
        MixedVec(i_offset + i) +=
            (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

        // scalar
        MixedVec(i_offset + i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau *
                                  scalar_sol[qp] * (*normals)[qp] * (*normals)[qp];

        // lm
        MixedVec(i_offset + i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * lm_sol[qp] *
                                  (*normals)[qp] * (*normals)[qp];
      }
    }
  }

  void scalar_face_jacobian(const unsigned int i_offset,
                            const unsigned int vector_j_offset,
                            const unsigned int scalar_j_offset,
                            const unsigned int lm_j_offset,
                            const unsigned int p_j_offset,
                            const unsigned int vel_component)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(i_offset + i, vector_j_offset + j) -=
              (*JxW_face)[qp] * mu * (*scalar_phi_face)[i][qp] *
              ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(p_n_dofs))
        {
          Gradient p_phi;
          p_phi(vel_component) = (*scalar_phi_face)[j][qp];
          // pressure
          MixedLM(i_offset + i, p_j_offset + j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
        }

        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(i_offset + i, scalar_j_offset + j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
              (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          MixedLM(i_offset + i, lm_j_offset + j) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] *
                                                    tau * (*lm_phi_face)[j][qp] * (*normals)[qp] *
                                                    (*normals)[qp];
      }
  }

  void lm_outlet_residual(const unsigned int i_offset,
                          const std::vector<Number> & scalar_sol,
                          const std::vector<Number> & lm_sol)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(lm_n_dofs))
      {
        // scalar
        LMVec(i_offset + i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * scalar_sol[qp] *
                               (*normals)[qp] * (*normals)[qp];

        // lm
        LMVec(i_offset + i) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * lm_sol[qp] *
                               (*normals)[qp] * (*normals)[qp];
      }
  }

  void lm_outlet_jacobian(const unsigned int i_offset,
                          const unsigned int scalar_j_offset,
                          const unsigned int lm_j_offset)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(lm_n_dofs))
      {
        for (const auto j : make_range(scalar_n_dofs))
          LMMixed(i_offset + i, scalar_j_offset + j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] *
                                                        tau * (*scalar_phi_face)[j][qp] *
                                                        (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          LMMat(i_offset + i, lm_j_offset + j) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau *
                                                  (*lm_phi_face)[j][qp] * (*normals)[qp] *
                                                  (*normals)[qp];
      }
  }

  void lm_face_residual(const unsigned int i_offset,
                        const std::vector<Gradient> & vector_sol,
                        const std::vector<Number> & scalar_sol,
                        const std::vector<Number> & lm_sol,
                        const unsigned int vel_component)
  {
    for (const auto qp : make_range(qface->n_points()))
    {
      Gradient qp_p;
      qp_p(vel_component) = p_sol[qp];
      for (const auto i : make_range(lm_n_dofs))
      {
        // vector
        LMVec(i_offset + i) -=
            (*JxW_face)[qp] * mu * (*lm_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // pressure
        LMVec(i_offset + i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

        // scalar
        LMVec(i_offset + i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * scalar_sol[qp] *
                               (*normals)[qp] * (*normals)[qp];

        // lm
        LMVec(i_offset + i) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * lm_sol[qp] *
                               (*normals)[qp] * (*normals)[qp];
      }
    }
  }

  void lm_face_jacobian(const unsigned int i_offset,
                        const unsigned int vector_j_offset,
                        const unsigned int scalar_j_offset,
                        const unsigned int lm_j_offset,
                        const unsigned int p_j_offset,
                        const unsigned int vel_component)
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(lm_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          LMMixed(i_offset + i, vector_j_offset + j) -=
              (*JxW_face)[qp] * mu * (*lm_phi_face)[i][qp] *
              ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(p_n_dofs))
        {
          Gradient p_phi;
          p_phi(vel_component) = (*scalar_phi_face)[j][qp];
          LMMat(i_offset + i, p_j_offset + j) +=
              (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
        }

        for (const auto j : make_range(scalar_n_dofs))
          LMMixed(i_offset + i, scalar_j_offset + j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] *
                                                        tau * (*scalar_phi_face)[j][qp] *
                                                        (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          LMMat(i_offset + i, lm_j_offset + j) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau *
                                                  (*lm_phi_face)[j][qp] * (*normals)[qp] *
                                                  (*normals)[qp];
      }
  }

  void assemble(const bool lm_solve)
  {
    auto & lm_soln_vector = lm_solve ? *lm_system->current_local_solution : *ghosted_old_solution;
    const auto u_num = mixed_system->variable_number("vel_x");
    const auto v_num = mixed_system->variable_number("vel_y");
    const auto qu_num = mixed_system->variable_number("qu");
    const auto qv_num = mixed_system->variable_number("qv");
    const auto lm_u_num = lm_system->variable_number("lm_u");
    const auto lm_v_num = lm_system->variable_number("lm_v");
    const auto p_num = lm_system->variable_number("pressure");

    std::vector<boundary_id_type> boundary_ids;
    const auto & boundary_info = mesh->get_boundary_info();

    for (const auto & elem : mesh->active_local_element_ptr_range())
    {
      // Retrive our dof indices for all fields
      mixed_dof_map->dof_indices(elem, qu_dof_indices, qu_num);
      mixed_dof_map->dof_indices(elem, u_dof_indices, u_num);
      lm_dof_map->dof_indices(elem, lm_u_dof_indices, lm_u_num);
      mixed_dof_map->dof_indices(elem, qv_dof_indices, qv_num);
      mixed_dof_map->dof_indices(elem, v_dof_indices, v_num);
      lm_dof_map->dof_indices(elem, lm_v_dof_indices, lm_v_num);
      lm_dof_map->dof_indices(elem, p_dof_indices, p_num);

      vector_n_dofs = qu_dof_indices.size();
      scalar_n_dofs = u_dof_indices.size();
      lm_n_dofs = lm_u_dof_indices.size();
      p_n_dofs = p_dof_indices.size();
      libmesh_assert(p_n_dofs == scalar_n_dofs);

      // Reinit our volume FE objects
      vector_fe->reinit(elem);
      scalar_fe->reinit(elem);

      libmesh_assert_equal_to(vector_n_dofs, vector_phi->size());
      libmesh_assert_equal_to(scalar_n_dofs, scalar_phi->size());

      const auto mixed_size = 2 * (vector_n_dofs + scalar_n_dofs);
      const auto lm_size = 2 * lm_n_dofs + p_n_dofs;

      // prepare our matrix/vector data structures
      MixedMat.setZero(mixed_size, mixed_size);
      MixedVec.setZero(mixed_size);
      LMMat.setZero(lm_size, lm_size);
      LMVec.setZero(lm_size);
      MixedLM.setZero(mixed_size, lm_size);
      LMMixed.setZero(lm_size, mixed_size);
      K_libmesh.resize(lm_size, lm_size);
      F_libmesh.resize(lm_size);

      // Get our local element dof values
      mixed_system->current_local_solution->get(qu_dof_indices, qu_dof_values);
      mixed_system->current_local_solution->get(u_dof_indices, u_dof_values);
      lm_soln_vector.get(lm_u_dof_indices, lm_u_dof_values);
      mixed_system->current_local_solution->get(qv_dof_indices, qv_dof_values);
      mixed_system->current_local_solution->get(v_dof_indices, v_dof_values);
      lm_soln_vector.get(lm_v_dof_indices, lm_v_dof_values);
      lm_soln_vector.get(p_dof_indices, p_dof_values);

      // Compute volumetric local qp solutions
      compute_qp_soln(qu_sol, qrule->n_points(), *vector_phi, qu_dof_values);
      compute_qp_soln(u_sol, qrule->n_points(), *scalar_phi, u_dof_values);
      compute_qp_soln(qv_sol, qrule->n_points(), *vector_phi, qv_dof_values);
      compute_qp_soln(v_sol, qrule->n_points(), *scalar_phi, v_dof_values);
      compute_qp_soln(p_sol, qrule->n_points(), *scalar_phi, p_dof_values);

      //
      // compute volumetric residuals and Jacobians
      //

      // qu and u
      vector_volume_residual(0, qu_sol, u_sol);
      scalar_volume_residual(vector_n_dofs, qu_sol, 0, sigma_u);
      vector_volume_jacobian(0, 0, vector_n_dofs);
      scalar_volume_jacobian(vector_n_dofs, 0, 2 * lm_n_dofs, 0);

      // qv and v
      vector_volume_residual(vector_n_dofs + scalar_n_dofs, qv_sol, v_sol);
      scalar_volume_residual(2 * vector_n_dofs + scalar_n_dofs, qv_sol, 1, sigma_v);
      vector_volume_jacobian(vector_n_dofs + scalar_n_dofs,
                             vector_n_dofs + scalar_n_dofs,
                             2 * vector_n_dofs + scalar_n_dofs);
      scalar_volume_jacobian(
          2 * vector_n_dofs + scalar_n_dofs, vector_n_dofs + scalar_n_dofs, 2 * lm_n_dofs, 1);

      // p
      pressure_volume_residual(2 * lm_n_dofs, u_sol, v_sol);
      pressure_volume_jacobian(2 * lm_n_dofs, vector_n_dofs, 2 * vector_n_dofs + scalar_n_dofs);

      for (auto side : elem->side_index_range())
      {
        // Reinit our face FE objects
        vector_fe_face->reinit(elem, side);
        scalar_fe_face->reinit(elem, side);
        lm_fe_face->reinit(elem, side);

        // Compute face local qp solutions
        compute_qp_soln(qu_sol, qface->n_points(), *vector_phi_face, qu_dof_values);
        compute_qp_soln(u_sol, qface->n_points(), *scalar_phi_face, u_dof_values);
        compute_qp_soln(lm_u_sol, qface->n_points(), *lm_phi_face, lm_u_dof_values);
        compute_qp_soln(qv_sol, qface->n_points(), *vector_phi_face, qv_dof_values);
        compute_qp_soln(v_sol, qface->n_points(), *scalar_phi_face, v_dof_values);
        compute_qp_soln(lm_v_sol, qface->n_points(), *lm_phi_face, lm_v_dof_values);
        compute_qp_soln(p_sol, qface->n_points(), *scalar_phi_face, p_dof_values);

        if (elem->neighbor_ptr(side) == nullptr)
        {
          boundary_info.boundary_ids(elem, side, boundary_ids);
          libmesh_assert(boundary_ids.size() == 1);
          const auto bnd_id = boundary_ids[0];
          if (bnd_id != right_bnd)
          {
            const auto dirichlet_velocity = [bnd_id, this]()
            {
              if (bnd_id == left_bnd)
                return RealVectorValue(1, 0);
              else
                return RealVectorValue(0, 0);
            }();

            // qu, u, lm_u
            vector_dirichlet_residual(0, dirichlet_velocity(0));
            scalar_dirichlet_residual(vector_n_dofs, qu_sol, u_sol, 0, dirichlet_velocity(0));
            scalar_dirichlet_jacobian(vector_n_dofs, 0, vector_n_dofs, 2 * lm_n_dofs, 0);

            // qv, v, lm_v
            vector_dirichlet_residual(vector_n_dofs + scalar_n_dofs, dirichlet_velocity(1));
            scalar_dirichlet_residual(
                2 * vector_n_dofs + scalar_n_dofs, qv_sol, v_sol, 1, dirichlet_velocity(1));
            scalar_dirichlet_jacobian(2 * vector_n_dofs + scalar_n_dofs,
                                      vector_n_dofs + scalar_n_dofs,
                                      2 * vector_n_dofs + scalar_n_dofs,
                                      2 * lm_n_dofs,
                                      1);

            // p
            pressure_dirichlet_residual(2 * lm_n_dofs, dirichlet_velocity);

            // Set the LMs on these Dirichlet boundary faces to 0
            create_identity_residual(*qface, *JxW_face, *lm_phi_face, lm_u_sol, lm_n_dofs, 0);
            create_identity_residual(
                *qface, *JxW_face, *lm_phi_face, lm_v_sol, lm_n_dofs, lm_n_dofs);
            create_identity_jacobian(*qface, *JxW_face, *lm_phi_face, lm_n_dofs, 0);
            create_identity_jacobian(*qface, *JxW_face, *lm_phi_face, lm_n_dofs, lm_n_dofs);

            continue;
          }
        }

        //
        // if we got here, then we are on an internal face or an outlet face
        //

        // qu, u, lm_u
        vector_face_residual(0, lm_u_sol);
        vector_face_jacobian(0, 0);
        scalar_face_residual(vector_n_dofs, qu_sol, u_sol, lm_u_sol, 0);
        scalar_face_jacobian(vector_n_dofs, 0, vector_n_dofs, 0, 2 * lm_n_dofs, 0);
        lm_face_residual(0, qu_sol, u_sol, lm_u_sol, 0);
        lm_face_jacobian(0, 0, vector_n_dofs, 0, 2 * lm_n_dofs, 0);

        // qv, v, lm_v
        vector_face_residual(vector_n_dofs + scalar_n_dofs, lm_v_sol);
        vector_face_jacobian(vector_n_dofs + scalar_n_dofs, lm_n_dofs);
        scalar_face_residual(2 * vector_n_dofs + scalar_n_dofs, qv_sol, v_sol, lm_v_sol, 1);
        scalar_face_jacobian(2 * vector_n_dofs + scalar_n_dofs,
                             vector_n_dofs + scalar_n_dofs,
                             2 * vector_n_dofs + scalar_n_dofs,
                             lm_n_dofs,
                             2 * lm_n_dofs,
                             1);
        lm_face_residual(lm_n_dofs, qv_sol, v_sol, lm_v_sol, 1);
        lm_face_jacobian(lm_n_dofs,
                         vector_n_dofs + scalar_n_dofs,
                         2 * vector_n_dofs + scalar_n_dofs,
                         lm_n_dofs,
                         2 * lm_n_dofs,
                         1);

        // p
        pressure_face_residual(2 * lm_n_dofs, lm_u_sol, lm_v_sol);
        pressure_face_jacobian(2 * lm_n_dofs, 0, lm_n_dofs);
      }

      MixedMatInv = MixedMat.inverse();
      lm_dof_indices = lm_u_dof_indices;
      lm_dof_indices.insert(lm_dof_indices.end(), lm_v_dof_indices.begin(), lm_v_dof_indices.end());
      lm_dof_indices.insert(lm_dof_indices.end(), p_dof_indices.begin(), p_dof_indices.end());
      libmesh_assert(lm_size == lm_dof_indices.size());

      if (lm_solve)
      {
        const auto LMProductMat = -LMMixed * MixedMatInv;
        LMMat += LMProductMat * MixedLM;
        LMVec += LMProductMat * MixedVec;
        libmesh_assert(cast_int<std::size_t>(LMMat.rows()) == lm_size);
        libmesh_assert(cast_int<std::size_t>(LMMat.cols()) == lm_size);
        libmesh_assert(cast_int<std::size_t>(LMVec.size()) == lm_size);

        for (const auto i : make_range(lm_size))
        {
          for (const auto j : make_range(lm_size))
            K_libmesh(i, j) = LMMat(i, j);
          F_libmesh(i) = LMVec(i);
        }

        // We were performing our finite element assembly for the implicit solve step of our
        // example. Add our local element vectors/matrices into the global system
        J->add_matrix(K_libmesh, lm_dof_indices);
        residual->add_vector(F_libmesh, lm_dof_indices);
      }
      else
      {
        //
        // We are doing our finite element assembly for the second time. We now know the Lagrange
        // multiplier solution. With that and the local element matrices and vectors we can compute
        // the vector and scalar solutions
        //
        ghosted_increment->get(lm_dof_indices, lm_increment_dof_values);
        LMIncrement.resize(lm_size);
        for (const auto i : index_range(lm_dof_indices))
          LMIncrement(i) = lm_increment_dof_values[i];

        MixedIncrement = MixedMatInv * (-MixedVec - MixedLM * LMIncrement);
        mixed_dof_indices = qu_dof_indices;
        auto append = [this](const auto & dofs)
        { mixed_dof_indices.insert(mixed_dof_indices.end(), dofs.begin(), dofs.end()); };
        append(u_dof_indices);
        append(qv_dof_indices);
        append(v_dof_indices);
        libmesh_assert(mixed_dof_indices.size() == mixed_size);
        mixed_increment_dof_values.resize(mixed_size);
        for (const auto i : index_range(mixed_increment_dof_values))
          mixed_increment_dof_values[i] = MixedIncrement(i);

        mixed_system->solution->add_vector(mixed_increment_dof_values, mixed_dof_indices);
      }
    }

    if (!lm_solve)
    {
      mixed_system->solution->close();
      // Scatter solution into the current_solution which is used in error computation
      mixed_system->update();
    }
  }

  // volume
  const std::vector<Real> * JxW;
  const std::vector<Point> * q_point;
  const std::vector<std::vector<RealVectorValue>> * vector_phi;
  const std::vector<std::vector<Real>> * scalar_phi;
  const std::vector<std::vector<RealVectorValue>> * grad_scalar_phi;
  const std::vector<std::vector<Real>> * div_vector_phi;

  // face
  const std::vector<Real> * JxW_face;
  const std::vector<Point> * qface_point;
  const std::vector<std::vector<RealVectorValue>> * vector_phi_face;
  const std::vector<std::vector<Real>> * scalar_phi_face;
  const std::vector<std::vector<Real>> * lm_phi_face;
  const std::vector<Point> * normals;

  //
  // We will need "Eigen" versions of many of the matrices/vectors because the
  // libMesh DenseMatrix doesn't have an inverse API. Let's match the notation of Nguyen's "An
  // implicit high-order hybridizable discontinuous Galerkin method for nonlinear
  // convection–diffusion equations" with the exception that instead of a fancy F we use P for
  // the residual vector for the LMs post-elimination
  //

  // LM matrix and RHS after eliminating vector and scalar dofs
  DenseMatrix<Number> K_libmesh;
  DenseVector<Number> F_libmesh;

  // On-diagonals
  EigenMatrix MixedMat, MixedMatInv, LMMat;
  // Vectors
  EigenVector MixedVec, LMVec;
  // Off-diagonals
  EigenMatrix MixedLM, LMMixed;

  // Containers for dof indices
  std::vector<dof_id_type> qu_dof_indices;
  std::vector<dof_id_type> u_dof_indices;
  std::vector<dof_id_type> lm_u_dof_indices;
  std::vector<dof_id_type> qv_dof_indices;
  std::vector<dof_id_type> v_dof_indices;
  std::vector<dof_id_type> lm_v_dof_indices;
  std::vector<dof_id_type> p_dof_indices;
  std::vector<dof_id_type> lm_dof_indices;
  std::vector<dof_id_type> mixed_dof_indices;

  // local solutions at quadrature points
  std::vector<Gradient> qu_sol;
  std::vector<Number> u_sol;
  std::vector<Number> lm_u_sol;
  std::vector<Gradient> qv_sol;
  std::vector<Number> v_sol;
  std::vector<Number> lm_v_sol;
  std::vector<Number> p_sol;

  // stresses at quadrature points
  std::vector<Gradient> sigma_u;
  std::vector<Gradient> sigma_v;

  // local degree of freedom values
  std::vector<Number> qu_dof_values;
  std::vector<Number> u_dof_values;
  std::vector<Number> lm_u_dof_values;
  std::vector<Number> qv_dof_values;
  std::vector<Number> v_dof_values;
  std::vector<Number> lm_v_dof_values;
  std::vector<Number> p_dof_values;

  // local degree of freedom increment values
  std::vector<Number> lm_increment_dof_values;
  std::vector<Number> mixed_increment_dof_values;
  EigenVector LMIncrement, MixedIncrement;

  // Number of dofs on elem
  std::size_t vector_n_dofs;
  std::size_t scalar_n_dofs;
  std::size_t lm_n_dofs;
  std::size_t p_n_dofs;

  // Our stabilization coefficient
  static constexpr Real tau = 1;

  // The viscosity
  static constexpr Real mu = 1;

  const USoln u_true_soln;
  const VSoln v_true_soln;
  const PSoln p_true_soln;
};

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Parse the input file.
  GetPot infile("vector_fe_ex9.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the command line and the input file.
  const unsigned int dimension = 2;
  const unsigned int grid_size = infile("grid_size", 2);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate first order side hierarchics, we must
  // use TRI6/7 elements
  const std::string elem_str = infile("element_type", std::string("TRI6"));

  libmesh_error_msg_if(elem_str != "TRI6" && elem_str != "TRI7",
                       "You selected "
                           << elem_str
                           << " but this example must be run with TRI6, TRI7, QUAD8, or QUAD9 in 2d"
                           << " or with TET14, or HEX27 in 3d.");

  MeshTools::Generation::build_square(
      mesh, 5 * grid_size, grid_size, 0., 10., -1, 1., Utility::string_to_enum<ElemType>(elem_str));

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the system "Mixed" and its variables.
  auto & system = equation_systems.add_system<System>("Mixed");

  // Add the LM system
  auto & lm_system = equation_systems.add_system<NonlinearImplicitSystem>("Lambda");

  // Adds the velocity variables and their gradients
  system.add_variable("qu", FIRST, L2_LAGRANGE_VEC);
  system.add_variable("qv", FIRST, L2_LAGRANGE_VEC);
  system.add_variable("vel_x", FIRST, L2_LAGRANGE);
  system.add_variable("vel_y", FIRST, L2_LAGRANGE);

  // Add our Lagrange multiplier to the implicit system
  lm_system.add_variable("lm_u", FIRST, SIDE_HIERARCHIC);
  lm_system.add_variable("lm_v", FIRST, SIDE_HIERARCHIC);
  lm_system.add_variable("pressure", FIRST, L2_LAGRANGE);

  // Add vectors for increment
  auto & ghosted_inc = lm_system.add_vector("ghosted_increment", true, GHOSTED);
  auto & parallel_inc = lm_system.add_vector("parallel_increment", true, PARALLEL);
  auto & ghosted_old = lm_system.add_vector("ghosted_old", true, GHOSTED);

  const FEType vector_fe_type(FIRST, L2_LAGRANGE_VEC);
  const FEType scalar_fe_type(FIRST, L2_LAGRANGE);
  const FEType lm_fe_type(FIRST, SIDE_HIERARCHIC);

  HDGProblem hdg;
  hdg.mesh = &mesh;
  hdg.mixed_system = &system;
  hdg.lm_system = &lm_system;
  hdg.mixed_dof_map = &system.get_dof_map();
  hdg.lm_dof_map = &lm_system.get_dof_map();
  hdg.vector_fe = FEVectorBase::build(dimension, vector_fe_type);
  hdg.scalar_fe = FEBase::build(dimension, scalar_fe_type);
  hdg.qrule = QBase::build(QGAUSS, dimension, FIFTH);
  hdg.qface = QBase::build(QGAUSS, dimension - 1, FIFTH);
  hdg.vector_fe_face = FEVectorBase::build(dimension, vector_fe_type);
  hdg.scalar_fe_face = FEBase::build(dimension, scalar_fe_type);
  hdg.lm_fe_face = FEBase::build(dimension, lm_fe_type);
  hdg.ghosted_increment = &ghosted_inc;
  hdg.parallel_increment = &parallel_inc;
  hdg.ghosted_old_solution = &ghosted_old;

  lm_system.nonlinear_solver->residual_and_jacobian_object = &hdg;
  lm_system.nonlinear_solver->postcheck_object = &hdg;

  hdg.init();

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Solve the implicit system for the Lagrange multiplier
  lm_system.solve();

  // //
  // // Now we will compute our solution approximation errors
  // //

  // ExactSolution exact_sol(equation_systems);

  // if (dimension == 2)
  // {
  //   SolutionFunction<2> soln_func;
  //   SolutionGradient<2> soln_grad;

  //   // Build FunctionBase* containers to attach to the ExactSolution object.
  //   std::vector<FunctionBase<Number> *> sols(1, &soln_func);
  //   std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

  //   exact_sol.attach_exact_values(sols);
  //   exact_sol.attach_exact_derivs(grads);
  // }
  // else if (dimension == 3)
  // {
  //   SolutionFunction<3> soln_func;
  //   SolutionGradient<3> soln_grad;

  //   // Build FunctionBase* containers to attach to the ExactSolution object.
  //   std::vector<FunctionBase<Number> *> sols(1, &soln_func);
  //   std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

  //   exact_sol.attach_exact_values(sols);
  //   exact_sol.attach_exact_derivs(grads);
  // }

  // // Use higher quadrature order for more accurate error results.
  // int extra_error_quadrature = infile("extra_error_quadrature", 2);
  // if (extra_error_quadrature)
  //   exact_sol.extra_quadrature_order(extra_error_quadrature);

  // // Compute the error.
  // exact_sol.compute_error("Mixed", "qu");
  // exact_sol.compute_error("Mixed", "u");
  // exact_sol.compute_error("Mixed", "qv");
  // exact_sol.compute_error("Mixed", "v");

  // // Print out the error values.
  // libMesh::out << "L2 error for qu is: " << exact_sol.l2_error("Mixed", "qu") << std::endl;
  // libMesh::out << "L2 error for u is: " << exact_sol.l2_error("Mixed", "u") << std::endl;
  // libMesh::out << "L2 error for qv is: " << exact_sol.l2_error("Mixed", "qv") << std::endl;
  // libMesh::out << "L2 error for v is: " << exact_sol.l2_error("Mixed", "v") << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}

#else

int
main()
{
  return 0;
}

#endif
#else

int
main()
{
  return 0;
}

#endif
