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

// Compute the stabilization parameter
Real
compute_tau(const bool internal_face, bool & tau_found, const Elem * const elem)
{
  if (!internal_face)
    // if we're an external face then tau is 0
    return 0;
  else if (tau_found)
    // if we've already applied a non-zero tau on another face, then we also return 0
    return 0;
  else
  {
    // This is our first internal face. We will apply our non-zero tau (and mark that we have)
    tau_found = true;
    return 1 / elem->hmin();
  }
}

class HDGProblem : public libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian,
                   public libMesh::NonlinearImplicitSystem::ComputePostCheck
{
public:
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
  void vector_volume_residual()
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(vector_n_dofs))
      {
        // Vector equation dependence on vector dofs
        R(i) -= (*JxW)[qp] * ((*vector_phi)[i][qp] * vector_sol[qp]);

        // Vector equation dependence on scalar dofs
        R(i) += (*JxW)[qp] * ((*div_vector_phi)[i][qp] * scalar_sol[qp]);
      }
  }

  void vector_volume_jacobian()
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(vector_n_dofs))
      {
        // Vector equation dependence on vector dofs
        for (const auto j : make_range(vector_n_dofs))
          A(i, j) -= (*JxW)[qp] * ((*vector_phi)[i][qp] * (*vector_phi)[j][qp]);

        // Vector equation dependence on scalar dofs
        for (const auto j : make_range(scalar_n_dofs))
          Bt(i, j) += (*JxW)[qp] * ((*div_vector_phi)[i][qp] * (*scalar_phi)[j][qp]);
      }
  }

  void scalar_volume_residual()
  {
    const auto dim = mesh->mesh_dimension();
    for (const auto qp : make_range(qrule->n_points()))
    {
      // Prepare forcing function
      const Real x = (*q_point)[qp](0);
      const Real y = (*q_point)[qp](1);
      const Real z = (*q_point)[qp](2);

      // "f" is the forcing function for the Poisson equation, which is
      // just the divergence of the exact solution for the vector field.
      // This is the well-known "method of manufactured solutions".
      Real f = 0;
      if (dim == 2)
        f = MixedExactSolution().forcing(x, y);
      else if (dim == 3)
        f = MixedExactSolution().forcing(x, y, z);

      for (const auto i : make_range(scalar_n_dofs))
      {
        // Scalar equation dependence on vector dofs
        F(i) += (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vector_sol[qp]);

        // Scalar equation RHS
        F(i) += (*JxW)[qp] * (*scalar_phi)[i][qp] * f;
      }
    }
  }

  void scalar_volume_jacobian()
  {
    for (const auto qp : make_range(qrule->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
        // Scalar equation dependence on vector dofs
        for (const auto j : make_range(vector_n_dofs))
          B(i, j) += (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * (*vector_phi)[j][qp]);
  }

  void vector_dirichlet_residual()
  {
    const auto dim = mesh->mesh_dimension();
    for (const auto qp : make_range(qface->n_points()))
    {
      const Real xf = (*qface_point)[qp](0);
      const Real yf = (*qface_point)[qp](1);
      const Real zf = (*qface_point)[qp](2);

      // The boundary value for scalar field.
      Real scalar_value = 0;
      if (dim == 2)
        scalar_value = MixedExactSolution().scalar(xf, yf);
      else if (dim == 3)
        scalar_value = MixedExactSolution().scalar(xf, yf, zf);

      // External boundary -> Dirichlet faces -> Vector equation RHS
      for (const auto i : make_range(vector_n_dofs))
        R(i) -= (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * scalar_value;
    }
  }

  void vector_face_residual()
  {
    for (const auto qp : make_range(qface->n_points()))
      // Vector equation dependence on LM dofs
      for (const auto i : make_range(vector_n_dofs))
        R(i) -= (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * lm_sol[qp];
  }

  void vector_face_jacobian()
  {
    for (const auto qp : make_range(qface->n_points()))
      // Vector equation dependence on LM dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(lm_n_dofs))
          Ct(i, j) -= (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) *
                      (*lm_phi_face)[j][qp];
  }

  void scalar_face_residual()
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        // vector
        F(i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // scalar
        F(i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * scalar_sol[qp] *
                (*normals)[qp] * (*normals)[qp];

        // lm
        F(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * lm_sol[qp] * (*normals)[qp] *
                (*normals)[qp];
      }
  }

  void scalar_face_jacobian()
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(scalar_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          B(i, j) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] *
                     ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(scalar_n_dofs))
          D(i, j) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
                     (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          E(i, j) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*lm_phi_face)[j][qp] *
                     (*normals)[qp] * (*normals)[qp];
      }
  }

  void lm_face_residual()
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(lm_n_dofs))
      {
        // vector
        L(i) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // scalar
        L(i) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * scalar_sol[qp] * (*normals)[qp] *
                (*normals)[qp];

        // lm
        L(i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * lm_sol[qp] * (*normals)[qp] *
                (*normals)[qp];
      }
  }

  void lm_face_jacobian()
  {
    for (const auto qp : make_range(qface->n_points()))
      for (const auto i : make_range(lm_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          C(i, j) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] *
                     ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(scalar_n_dofs))
          G(i, j) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
                     (*normals)[qp] * (*normals)[qp];

        for (const auto j : make_range(lm_n_dofs))
          H(i, j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * (*lm_phi_face)[j][qp] *
                     (*normals)[qp] * (*normals)[qp];
      }
  }

  void assemble(const bool lm_solve)
  {
    auto & lm_soln_vector = lm_solve ? *lm_system->current_local_solution : *ghosted_old_solution;
    for (const auto & elem : mesh->active_local_element_ptr_range())
    {
      // Retrive our dof indices for all fields
      mixed_dof_map->dof_indices(elem, vector_dof_indices, mixed_system->variable_number("q"));
      mixed_dof_map->dof_indices(elem, scalar_dof_indices, mixed_system->variable_number("u"));
      lm_dof_map->dof_indices(elem, lm_dof_indices, lm_system->variable_number("lambda"));

      vector_n_dofs = vector_dof_indices.size();
      scalar_n_dofs = scalar_dof_indices.size();
      lm_n_dofs = lm_dof_indices.size();

      // Reinit our volume FE objects
      vector_fe->reinit(elem);
      scalar_fe->reinit(elem);

      libmesh_assert_equal_to(vector_n_dofs, vector_phi->size());
      libmesh_assert_equal_to(scalar_n_dofs, scalar_phi->size());

      // prepare our matrix/vector data structures for the vector equation
      A.setZero(vector_n_dofs, vector_n_dofs);
      Bt.setZero(vector_n_dofs, scalar_n_dofs);
      Ct.setZero(vector_n_dofs, lm_n_dofs);
      R.setZero(vector_n_dofs);
      // and for the scalar equation
      B.setZero(scalar_n_dofs, vector_n_dofs);
      D.setZero(scalar_n_dofs, scalar_n_dofs);
      E.setZero(scalar_n_dofs, lm_n_dofs);
      F.setZero(scalar_n_dofs);
      // and for the LM equation
      C.setZero(lm_n_dofs, vector_n_dofs);
      G.setZero(lm_n_dofs, scalar_n_dofs);
      H.setZero(lm_n_dofs, lm_n_dofs);
      L.setZero(lm_n_dofs);

      // Get our local element dof values
      mixed_system->current_local_solution->get(vector_dof_indices, vector_dof_values);
      mixed_system->current_local_solution->get(scalar_dof_indices, scalar_dof_values);
      lm_soln_vector.get(lm_dof_indices, lm_dof_values);

      // Compute volumetric local qp solutions
      compute_qp_soln(vector_sol, qrule->n_points(), *vector_phi, vector_dof_values);
      compute_qp_soln(scalar_sol, qrule->n_points(), *scalar_phi, scalar_dof_values);

      // compute volumetric residuals and Jacobians
      vector_volume_residual();
      scalar_volume_residual();
      vector_volume_jacobian();
      scalar_volume_jacobian();

      // At the beginning of the loop, we mark that we haven't found our "Single-Face" yet
      bool tau_found = false;
      for (auto side : elem->side_index_range())
      {
        // Reinit our face FE objects
        vector_fe_face->reinit(elem, side);
        scalar_fe_face->reinit(elem, side);
        lm_fe_face->reinit(elem, side);

        // Compute face local qp solutions
        compute_qp_soln(vector_sol, qface->n_points(), *vector_phi_face, vector_dof_values);
        compute_qp_soln(scalar_sol, qface->n_points(), *scalar_phi_face, scalar_dof_values);
        compute_qp_soln(lm_sol, qface->n_points(), *lm_phi_face, lm_dof_values);

        tau = compute_tau(elem->neighbor_ptr(side) != nullptr, tau_found, elem);

        if (elem->neighbor_ptr(side) == nullptr)
        {
          vector_dirichlet_residual();
          scalar_face_residual();
          scalar_face_jacobian();
          for (const auto qp : make_range(qface->n_points()))
          {
            // Need to do something with the external boundary LM dofs to prevent
            // the matrix from being singular
            for (const auto i : make_range(lm_n_dofs))
            {
              L(i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * lm_sol[qp];
              for (const auto j : make_range(lm_n_dofs))
                H(i, j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (*lm_phi_face)[j][qp];
            }
          }
        }
        else
        {
          vector_face_residual();
          vector_face_jacobian();
          scalar_face_residual();
          scalar_face_jacobian();
          lm_face_residual();
          lm_face_jacobian();
        }
      }

      const auto mixed_size = vector_n_dofs + scalar_n_dofs;
      MixedMat.resize(mixed_size, mixed_size);
      MixedLM.resize(mixed_size, lm_n_dofs);
      MixedVec.resize(mixed_size);
      for (const auto i : make_range(vector_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(i, j) = A(i, j);

        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(i, vector_n_dofs + j) = Bt(i, j);

        for (const auto j : make_range(lm_n_dofs))
          MixedLM(i, j) = Ct(i, j);

        MixedVec(i) = R(i);
      }
      for (const auto i : make_range(scalar_n_dofs))
      {
        for (const auto j : make_range(vector_n_dofs))
          MixedMat(vector_n_dofs + i, j) = B(i, j);

        for (const auto j : make_range(scalar_n_dofs))
          MixedMat(vector_n_dofs + i, vector_n_dofs + j) = D(i, j);

        for (const auto j : make_range(lm_n_dofs))
          MixedLM(vector_n_dofs + i, j) = E(i, j);

        MixedVec(vector_n_dofs + i) = F(i);
      }

      MixedMatInv = MixedMat.inverse();

      if (lm_solve)
      {
        LMMixed.resize(lm_n_dofs, mixed_size);
        for (const auto i : make_range(lm_n_dofs))
        {
          for (const auto j : make_range(vector_n_dofs))
            LMMixed(i, j) = C(i, j);

          for (const auto j : make_range(scalar_n_dofs))
            LMMixed(i, vector_n_dofs + j) = G(i, j);
        }

        const auto LMProductMat = -LMMixed * MixedMatInv;
        K = LMProductMat * MixedLM + H;
        P = LMProductMat * MixedVec + L;
        K_libmesh.resize(lm_n_dofs, lm_n_dofs);
        P_libmesh.resize(lm_n_dofs);

        for (const auto i : make_range(lm_n_dofs))
        {
          for (const auto j : make_range(lm_n_dofs))
            K_libmesh(i, j) = K(i, j);
          P_libmesh(i) = P(i);
        }

        // We were performing our finite element assembly for the implicit solve step of our
        // example. Add our local element vectors/matrices into the global system
        J->add_matrix(K_libmesh, lm_dof_indices);
        residual->add_vector(P_libmesh, lm_dof_indices);
      }
      else
      {
        //
        // We are doing our finite element assembly for the second time. We now know the Lagrange
        // multiplier solution. With that and the local element matrices and vectors we can compute
        // the vector and scalar solutions
        //
        ghosted_increment->get(lm_dof_indices, lm_increment_dof_values);
        LMIncrement.resize(lm_n_dofs);
        for (const auto i : make_range(lm_n_dofs))
          LMIncrement(i) = lm_increment_dof_values[i];

        MixedIncrement = MixedMatInv * (-MixedVec - MixedLM * LMIncrement);
        vector_increment_dof_values.resize(vector_n_dofs);
        scalar_increment_dof_values.resize(scalar_n_dofs);
        for (const auto i : make_range(vector_n_dofs))
          vector_increment_dof_values[i] = MixedIncrement(i);
        for (const auto i : make_range(scalar_n_dofs))
          scalar_increment_dof_values[i] = MixedIncrement(vector_n_dofs + i);

        mixed_system->solution->add_vector(vector_increment_dof_values, vector_dof_indices);
        mixed_system->solution->add_vector(scalar_increment_dof_values, scalar_dof_indices);
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
  DenseVector<Number> P_libmesh;
  EigenMatrix K;
  EigenVector P;
  // (A  Bt Ct)(q) = (R)
  // (B  D  E) (u) = (F)
  // (C  G  H) (l) = (L)
  EigenMatrix A, Bt, Ct, B, D, E, C, G, H;
  EigenVector R, F, L;

  // Compositions
  EigenMatrix MixedMat, MixedMatInv;
  EigenVector MixedVec;
  EigenMatrix MixedLM, LMMixed;

  // Containers for dof indices
  std::vector<dof_id_type> vector_dof_indices;
  std::vector<dof_id_type> scalar_dof_indices;
  std::vector<dof_id_type> lm_dof_indices;

  // local solutions at quadrature points
  std::vector<Gradient> vector_sol;
  std::vector<Number> scalar_sol;
  std::vector<Number> lm_sol;

  // local degree of freedom values
  std::vector<Number> vector_dof_values;
  std::vector<Number> scalar_dof_values;
  std::vector<Number> lm_dof_values;

  // local degree of freedom increment values
  std::vector<Number> vector_increment_dof_values;
  std::vector<Number> scalar_increment_dof_values;
  std::vector<Number> lm_increment_dof_values;
  EigenVector LMIncrement, MixedIncrement;

  // Number of dofs on elem
  std::size_t vector_n_dofs;
  std::size_t scalar_n_dofs;
  std::size_t lm_n_dofs;

  // Our stabilization coefficient
  Real tau;
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
  const unsigned int dimension = infile("dim", 2);
  const unsigned int grid_size = infile("grid_size", 15);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate first order side hierarchics, we must
  // use TRI6/7 or QUAD8/9 elements in 2d, or TET14 or HEX27 in 3d.
  const std::string elem_str = infile("element_type", std::string("TRI6"));

  libmesh_error_msg_if((dimension == 2 && elem_str != "TRI6" && elem_str != "TRI7" &&
                        elem_str != "QUAD8" && elem_str != "QUAD9") ||
                           (dimension == 3 && elem_str != "TET14" && elem_str != "HEX27"),
                       "You selected "
                           << elem_str
                           << " but this example must be run with TRI6, TRI7, QUAD8, or QUAD9 in 2d"
                           << " or with TET14, or HEX27 in 3d.");

  if (dimension == 2)
    MeshTools::Generation::build_square(
        mesh, grid_size, grid_size, -1., 1., -1., 1., Utility::string_to_enum<ElemType>(elem_str));
  else if (dimension == 3)
    MeshTools::Generation::build_cube(mesh,
                                      grid_size,
                                      grid_size,
                                      grid_size,
                                      -1.,
                                      1.,
                                      -1.,
                                      1.,
                                      -1.,
                                      1.,
                                      Utility::string_to_enum<ElemType>(elem_str));

  // Make sure the code is robust against nodal reorderings.
  MeshTools::Modification::permute_elements(mesh);

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the system "Mixed" and its variables.
  auto & system = equation_systems.add_system<System>("Mixed");

  // Add the LM system
  auto & lm_system = equation_systems.add_system<NonlinearImplicitSystem>("Lambda");

  // Adds the variable "q" and "u" to "Mixed". "q" will be our vector field
  // whereas "u" will be the scalar field.
  system.add_variable("q", FIRST, L2_LAGRANGE_VEC);
  system.add_variable("u", FIRST, L2_LAGRANGE);

  // Add our Lagrange multiplier to the implicit system
  lm_system.add_variable("lambda", FIRST, SIDE_HIERARCHIC);

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

  //
  // Now we will compute our solution approximation errors
  //

  ExactSolution exact_sol(equation_systems);

  if (dimension == 2)
  {
    SolutionFunction<2> soln_func;
    SolutionGradient<2> soln_grad;

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }
  else if (dimension == 3)
  {
    SolutionFunction<3> soln_func;
    SolutionGradient<3> soln_grad;

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }

  // Use higher quadrature order for more accurate error results.
  int extra_error_quadrature = infile("extra_error_quadrature", 2);
  if (extra_error_quadrature)
    exact_sol.extra_quadrature_order(extra_error_quadrature);

  // Compute the error.
  exact_sol.compute_error("Mixed", "q");
  exact_sol.compute_error("Mixed", "u");

  // Print out the error values.
  libMesh::out << "L2 error is: " << exact_sol.l2_error("Mixed", "q") << std::endl;
  libMesh::out << "L2 error for u is: " << exact_sol.l2_error("Mixed", "u") << std::endl;

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
