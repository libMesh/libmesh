// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/variational_smoother_system.h"

#include "libmesh/elem.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"
#include "libmesh/enum_to_string.h"
#include <libmesh/reference_elem.h>

// C++ includes
#include <functional> // std::reference_wrapper

namespace libMesh
{

/**
 * Function to prevent dividing by zero for degenerate elements
 */
Real chi_epsilon(const Real & x, const Real epsilon_squared)
{
  return 0.5 * (x + std::sqrt(epsilon_squared + Utility::pow<2>(x)));
}

/**
 * Given an fe_map, element dimension, and quadrature point index, returns the
 * Jacobian of the physical-to-reference mapping.
 */
RealTensor get_jacobian_at_qp(const FEMap & fe_map,
                          const unsigned int & dim,
                          const unsigned int & qp)
{
  const auto & dxyzdxi = fe_map.get_dxyzdxi()[qp];
  const auto & dxyzdeta = fe_map.get_dxyzdeta()[qp];
  const auto & dxyzdzeta = fe_map.get_dxyzdzeta()[qp];

  // RealTensors are always 3x3, so we will fill any dimensions above dim
  // with 1s on the diagonal. This indicates a 1 to 1 relationship between
  // the physical and reference elements in these extra dimensions.
  switch (dim)
  {
    case 1:
      return RealTensor(
        dxyzdxi(0), 0, 0,
                     0, 1, 0,
                     0, 0, 1
      );
      break;

    case 2:
      return RealTensor(
        dxyzdxi(0), dxyzdeta(0), 0,
        dxyzdxi(1), dxyzdeta(1), 0,
        0,              0,               1
      );
      break;

    case 3:
      return RealTensor(
        dxyzdxi,
        dxyzdeta,
        dxyzdzeta
      ).transpose();  // Note the transposition!
      break;

    default:
      libmesh_error_msg("Unsupported dimension.");
  }
}

/**
 * Compute the trace of a dim-dimensional matrix.
 */
Real trace(const RealTensor & A, const unsigned int & dim)
{
  Real tr = 0.0;
  for (const auto i : make_range(dim))
    tr += A(i, i);

  return tr;
}

VariationalSmootherSystem::~VariationalSmootherSystem () = default;

void VariationalSmootherSystem::assembly (bool get_residual,
                                          bool get_jacobian,
                                          bool apply_heterogeneous_constraints,
                                          bool apply_no_constraints)
{
  // Update the mesh based on the current variable values
  auto & mesh = this->get_mesh();
  this->solution->close();

  for (auto * node : mesh.local_node_ptr_range())
  {
    for (const auto d : make_range(mesh.mesh_dimension()))
    {
      const auto dof_id = node->dof_number(this->number(), d, 0);
      // Update mesh
      (*node)(d) = libmesh_real((*current_local_solution)(dof_id));
    }
  }

  SyncNodalPositions sync_object(mesh);
  Parallel::sync_dofobject_data_by_id (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync_object);

  // Compute and update mesh quality information
  compute_mesh_quality_info();

  // Update _epsilon_squared_assembly based on whether the mesh is tangled
  _epsilon_squared_assembly = _mesh_info.mesh_is_tangled ? _epsilon_squared : 0.;

  FEMSystem::assembly(get_residual, get_jacobian, apply_heterogeneous_constraints, apply_no_constraints);
}

void VariationalSmootherSystem::init_data ()
{
  auto & mesh = this->get_mesh();
  const auto & elem_orders = mesh.elem_default_orders();
  libmesh_error_msg_if(elem_orders.size() != 1,
      "The variational smoother cannot be used for mixed-order meshes!");
  const auto fe_order = *elem_orders.begin();
  // Add a variable for each dimension of the mesh
  // "r0" for x, "r1" for y, "r2" for z
  for (const auto & d : make_range(mesh.mesh_dimension()))
    this->add_variable ("r" + std::to_string(d), fe_order);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Set the current_local_solution to the current mesh for the initial guess
  this->solution->close();

  for (auto * node : mesh.local_node_ptr_range())
  {
    for (const auto d : make_range(mesh.mesh_dimension()))
    {
      const auto dof_id = node->dof_number(this->number(), d, 0);
      // Update solution
      (*solution).set(dof_id, (*node)(d));
    }
  }

  this->prepare_for_smoothing();
}

void VariationalSmootherSystem::prepare_for_smoothing()
{
  // If this method has already been called to set _ref_vol, just return.
  if (std::abs(_ref_vol) > TOLERANCE * TOLERANCE)
    return;

  std::unique_ptr<DiffContext> con = this->build_context();
  FEMContext & femcontext = cast_ref<FEMContext &>(*con);
  this->init_context(femcontext);

  const auto & mesh = this->get_mesh();

  Real elem_averaged_det_S_sum = 0.;

  // Make pre-requests before reinit() for efficiency in
  // --enable-deprecated builds, and to avoid errors in
  // --disable-deprecated builds.
  const auto & fe_map = femcontext.get_element_fe(0)->get_fe_map();
  const auto & JxW = fe_map.get_JxW();

  for (const auto * elem : mesh.active_local_element_ptr_range())
    {
      femcontext.pre_fe_reinit(*this, elem);
      femcontext.elem_fe_reinit();

      // Add target element info, if applicable
      if (_target_jacobians.find(elem->type()) == _target_jacobians.end())
        {
          const auto [target_elem, target_nodes] = get_target_elem(elem->type());
          get_target_to_reference_jacobian(target_elem.get(),
                                           femcontext,
                                           _target_jacobians[elem->type()],
                                           _target_jacobian_dets[elem->type()]);
        }// if find == end()

      // Reference volume computation
      Real elem_integrated_det_S = 0.;
      for (const auto qp : index_range(JxW))
        elem_integrated_det_S += JxW[qp] * _target_jacobian_dets[elem->type()][qp];
      const auto ref_elem_vol = elem->reference_elem()->volume();
      elem_averaged_det_S_sum += elem_integrated_det_S / ref_elem_vol;

    } // for elem

  // Get contributions from elements on other processors
  mesh.comm().sum(elem_averaged_det_S_sum);

  _ref_vol = elem_averaged_det_S_sum / mesh.n_active_elem();
}

void VariationalSmootherSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * my_fe = nullptr;

  // Now make sure we have requested all the data
  // we need to build the system.

  // We might have a multi-dimensional mesh
  const std::set<unsigned char> & elem_dims =
    c.elem_dimensions();

  for (const auto & dim : elem_dims)
    {
      c.get_element_fe( 0, my_fe, dim );
      my_fe->get_nothing();

      auto & fe_map = my_fe->get_fe_map();
      fe_map.get_dxyzdxi();
      fe_map.get_dxyzdeta();
      fe_map.get_dxyzdzeta();
      fe_map.get_JxW();

      // Mesh may be tangled, allow negative Jacobians
      fe_map.set_jacobian_tolerance(std::numeric_limits<Real>::lowest());

      c.get_side_fe( 0, my_fe, dim );
      my_fe->get_nothing();
    }

  FEMSystem::init_context(context);
}


bool VariationalSmootherSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const Elem & elem = c.get_elem();

  unsigned int dim = c.get_dim();

  unsigned int x_var = 0, y_var = 1, z_var = 2;
  // In the case of lower dimensions, we will not access z (1D/2D) or y (1D)
  if (dim < 3)
  {
    z_var = 0;
    if (dim < 2)
      y_var = 0;
  }

  // The subvectors and submatrices we need to fill:
  // system residual
  std::reference_wrapper<DenseSubVector<Number>> F[3] =
  {
    c.get_elem_residual(x_var),
    c.get_elem_residual(y_var),
    c.get_elem_residual(z_var)
  };
  // system jacobian
  std::reference_wrapper<DenseSubMatrix<Number>> K[3][3] =
    {
      {c.get_elem_jacobian(x_var, x_var), c.get_elem_jacobian(x_var, y_var), c.get_elem_jacobian(x_var, z_var)},
      {c.get_elem_jacobian(y_var, x_var), c.get_elem_jacobian(y_var, y_var), c.get_elem_jacobian(y_var, z_var)},
      {c.get_elem_jacobian(z_var, x_var), c.get_elem_jacobian(z_var, y_var), c.get_elem_jacobian(z_var, z_var)}
    };

  // Quadrature info
  const auto quad_weights = c.get_element_qrule().get_weights();

  const auto distortion_weight = 1. - _dilation_weight;

  // Get some references to cell-specific data that
  // will be used to assemble the linear system.

  const auto & fe_map = c.get_element_fe(0)->get_fe_map();

  const auto & dphidxi_map = fe_map.get_dphidxi_map();
  const auto & dphideta_map = fe_map.get_dphideta_map();
  const auto & dphidzeta_map = fe_map.get_dphidzeta_map();

  // Integrate the distortion-dilation metric over the reference element
  for (const auto qp : index_range(quad_weights))
    {
      // Compute quantities needed to evaluate the distortion-dilation metric
      // and its gradient and Hessian.
      // The metric will be minimized when it's gradient with respect to the node
      // locations (R) is zero. For Newton's method, minimizing the gradient
      // requires computation of the gradient's Jacobian with respect to R.
      // The Jacobian of the distortion-dilation metric's gradientis the Hessian
      // of the metric.
      //
      // Note that the term "Jacobian" has two meanings in this mesh smoothing
      // application. The first meaning refers to the Jacobian w.r.t R of the
      // gradient w.r.t. R, (i.e., the Hessian of the metric). This is the K
      // variable defined above. This is also the Jacobian the 'request_jacobian'
      // variable refers to. The second mesning refers to the Jacobian of the
      // physical-to-reference element mapping. This is the Jacobian used to
      // compute the distorion-dilation metric.
      //
      // Grab the physical-to-reference mapping Jacobian matrix (i.e., "S") at this qp
      RealTensor S = get_jacobian_at_qp(fe_map, dim, qp);

      // Apply target element transformation
      S *= _target_jacobians[elem.type()][qp];

      // Compute quantities needed for the smoothing algorithm

      // determinant
      const Real det = S.det();
      const Real det_sq = det * det;
      const Real det_cube = det_sq * det;

      const Real ref_vol_sq = _ref_vol * _ref_vol;

      // trace of S^T * S
      // DO NOT USE RealTensor.tr for the trace, it will NOT be correct for
      // 1D and 2D meshes because of our hack of putting 1s in the diagonal of
      // S for the extra dimensions (see get_jacobian_at_qp)
      const auto tr = trace(S.transpose() * S, dim);
      const Real tr_div_dim = tr / dim;

      // Precompute pow(tr_div_dim, 0.5 * dim - x) for x = 0, 1, 2
      const Real half_dim = 0.5 * dim;
      const std::vector<Real> trace_powers{
        std::pow(tr_div_dim, half_dim),
        std::pow(tr_div_dim, half_dim - 1.),
        std::pow(tr_div_dim, half_dim - 2.),
      };

      // inverse of S
      const RealTensor S_inv = S.inverse();
      // inverse transpose of S
      const RealTensor S_inv_T = S_inv.transpose();

      // Identity matrix
      const RealTensor I(1, 0, 0,
                         0, 1, 0,
                         0, 0, 1);

      // The chi function allows us to handle degenerate elements
      const auto chi = chi_epsilon(det, _epsilon_squared_assembly);
      const Real chi_sq = chi * chi;
      const Real sqrt_term = std::sqrt(_epsilon_squared_assembly + det_sq);
      // dchi(x) / dx
      const Real chi_prime = 0.5 * (1. + det / sqrt_term);
      const Real chi_prime_sq = chi_prime * chi_prime;
      // d2chi(x) / dx2
      const Real chi_2prime = 0.5 * (1. / sqrt_term - det_sq / Utility::pow<3>(sqrt_term));

      // Distortion metric (beta)
      //const Real beta = trace_powers[0] / chi;
      const RealTensor dbeta_dS = (trace_powers[1] / chi) * S - (trace_powers[0] / chi_sq * chi_prime * det) * S_inv_T;

      // Dilation metric (mu)
      //const Real mu = 0.5 * (_ref_vol + det_sq / _ref_vol) / chi;
      // We represent d mu / dS as alpha(S) * S^-T, where alpha is a scalar function
      const Real alpha = (-chi_prime * det_cube + 2. * det_sq * chi - ref_vol_sq * det * chi_prime) / (2. * _ref_vol * chi_sq);
      const RealTensor dmu_dS = alpha * S_inv_T;

      // Combined metric (E)
      //const Real E = distortion_weight * beta + _dilation_weight * mu;
      const RealTensor dE_dS = distortion_weight * dbeta_dS + _dilation_weight * dmu_dS;

      // This vector is useful in computing dS/dR below
      std::vector<std::vector<std::vector<Real>>> dphi_maps = {dphidxi_map, dphideta_map, dphidzeta_map};

      // Compute residual (i.e., the gradient of the combined metric w.r.t node locations)
      // Recall that when the gradient (residual) is zero, the combined metric is minimized
      for (const auto l : elem.node_index_range())
      {
        for (const auto var_id : make_range(dim))
        {
          // Build dS/dR, the derivative of the physical-to-reference mapping Jacobin w.r.t. the mesh node locations
          RealTensor dS_dR = RealTensor(0);
          for (const auto jj : make_range(dim))
            dS_dR(var_id, jj) = dphi_maps[jj][l][qp];

          // Residual contribution. The contraction of dE/dS and dS/dR gives us
          // the gradient we are looking for, dE/dR
          F[var_id](l) += quad_weights[qp] * dE_dS.contract(dS_dR);
        }// for var_id
      }// for l


      if (request_jacobian)
      {
        // Compute jacobian of the smoothing system (i.e., the Hessian of the
        // combined metric w.r.t. the mesh node locations)

        // Precompute coefficients to be applied to each component of tensor
        // products in the loops below. At first glance, these coefficients look
        // like gibberish, but everything in the Hessian has been verified by
        // taking finite differences of the gradient. We should probably write
        // down the derivations of the gradient and Hessian somewhere...

        // Recall that above, dbeta_dS takes the form:
        // d(beta)/dS = c1(S) * S - c2(S) * S_inv_T,
        // where c1 and c2 are scalar-valued functions.
        const std::vector<Real> d2beta_dS2_coefs_times_distortion_weight = {
          //Part 1: scaler coefficients of d(c1 * S) / dS
          //
          // multiplies I[i,a] x I[j,b]
          (trace_powers[1] / chi) * distortion_weight,
          // multiplies S[a,b] x S[i,j]
          (((dim - 2.) / dim) * trace_powers[2] / chi) * distortion_weight,
          // multiplies S_inv[b,a] * S[i,j]
          (-(trace_powers[1] / chi_sq) * chi_prime * det) * distortion_weight,
          //
          //Part 2: scaler coefficients of d(-c2 * S_inv_T) / dS
          //
          // multiplies S[a,b] x S_inv[j,i]
          (-(trace_powers[1] / chi_sq) * chi_prime * det) * distortion_weight,
          // multiplies S_inv[b,a] x S_inv[j,i]
          (trace_powers[0] * (det / chi_sq)
            * ((2. * chi_prime_sq / chi - chi_2prime) * det - chi_prime)) * distortion_weight,
          // multiplies S_inv[b,i] x S_inv[j,a]
          ((trace_powers[0] / chi_sq) * chi_prime * det) * distortion_weight,
        };

        // d alpha / dS has the form c(S) * S^-T, where c is the scalar coefficient defined below
        const Real dalpha_dS_coef_times_dilation_weight = ((det / (2. * _ref_vol * chi_sq))
          * (-4. * _ref_vol * alpha * chi * chi_prime - chi_2prime * det_cube
             - chi_prime * det_sq + 4 * chi * det - ref_vol_sq * (chi_prime + det * chi_2prime))) * _dilation_weight;

        // This is also useful to precompute
        const Real alpha_times_dilation_weight = alpha * _dilation_weight;

        /*

        To increase the efficiency of the Jacobian computation, we take
        advantage of l-p symmetry, ij-ab symmetry, and the sparsity pattern of
        dS_dR. We also factor all possible multipliers out of inner loops for
        efficiency. The result is code that is more difficult to read. For
        clarity, consult the pseudo-code below in this comment.

        for (const auto l: elem.node_index_range()) // Contribution to Hessian
        from node l
        {
          for (const auto var_id1 : make_range(dim)) // Contribution from each
        x/y/z component of node l
          {
            // Build dS/dR_l, the derivative of the physical-to-reference
        mapping Jacobin w.r.t.
            // the l-th node
            RealTensor dS_dR_l = RealTensor(0);
            for (const auto ii : make_range(dim))
              dS_dR_l(var_id1, ii) = dphi_maps[ii][l][qp];

            for (const auto p: elem.node_index_range()) // Contribution to
        Hessian from node p
            {
              for (const auto var_id2 : make_range(dim)) // Contribution from
        each x/y/z component of node p
              {
                // Build dS/dR_l, the derivative of the physical-to-reference
        mapping Jacobin w.r.t.
                // the p-th node
                RealTensor dS_dR_p = RealTensor(0);
                for (const auto jj : make_range(dim))
                  dS_dR_p(var_id2, jj) = dphi_maps[jj][p][qp];

                Real d2beta_dR2 = 0.;
                Real d2mu_dR2 = 0.;
                // Perform tensor contraction
                for (const auto i : make_range(dim))
                {
                  for (const auto j : make_range(dim))
                  {
                    for (const auto a : make_range(dim))
                    {
                      for (const auto b : make_range(dim))
                      {
                        // Nasty tensor products to be multiplied by
        d2beta_dS2_coefs to get d2(beta) / dS2 const std::vector<Real>
        d2beta_dS2_tensor_contributions =
                        {
                          I(i,a)     * I(j,b),
                          S(a,b)     * S(i,j),
                          S_inv(b,a) * S(i,j),
                          S(a,b)     * S_inv(j,i),
                          S_inv(b,a) * S_inv(j,i),
                          S_inv(j,a) * S_inv(b,i),
                        };

                        // Combine precomputed coefficients with tensor products
        to get d2(beta) / dS2 Real d2beta_dS2 = 0.; for (const auto comp_id :
        index_range(d2beta_dS2_coefs))
                        {
                          const Real contribution = d2beta_dS2_coefs[comp_id] *
        d2beta_dS2_tensor_contributions[comp_id]; d2beta_dS2 += contribution;
                        }

                        // Incorporate tensor product portion to get d2(mu) /
        dS2 const Real d2mu_dS2 = dalpha_dS_coef * S_inv(b,a) * S_inv(j,i) -
        alpha * S_inv(b,i) * S_inv(j,a);

                        // Chain rule to change d/dS to d/dR
                        d2beta_dR2 += d2beta_dS2 * dS_dR_l(a, b) * dS_dR_p(i,
        j); d2mu_dR2 += d2mu_dS2 * dS_dR_l(a, b) * dS_dR_p(i, j);

                      }// for b
                    }// for a
                  }// for j
                }// for i, end tensor contraction

                // Jacobian contribution
                K[var_id1][var_id2](l, p) += quad_weights[qp] * ((1. -
        _dilation_weight) * d2beta_dR2 + _dilation_weight * d2mu_dR2);

              }// for var_id2
            }// for p
          }// for var_id1
        }// for l

        End pseudo-code, begin efficient code

        */

        for (const auto l: elem.node_index_range()) // Contribution to Hessian from node l
        {
          for (const auto var_id1 : make_range(dim)) // Contribution from each x/y/z component of node l
          {
            // Build dS/dR_l, the derivative of the physical-to-reference mapping Jacobin w.r.t.
            // the l-th node
            RealTensor dS_dR_l = RealTensor(0);
            for (const auto ii : make_range(dim))
              dS_dR_l(var_id1, ii) = dphi_maps[ii][l][qp];

            // Jacobian is symmetric, only need to loop over lower triangular portion
            for (const auto p: make_range(l + 1)) // Contribution to Hessian from node p
            {
              for (const auto var_id2 : make_range(dim)) // Contribution from each x/y/z component of node p
              {
                // Build dS/dR_l, the derivative of the physical-to-reference mapping Jacobin w.r.t.
                // the p-th node
                RealTensor dS_dR_p = RealTensor(0);
                for (const auto jj : make_range(dim))
                  dS_dR_p(var_id2, jj) = dphi_maps[jj][p][qp];

                Real d2E_dR2 = 0.;
                // Perform tensor contraction
                for (const auto i : make_range(dim))
                {

                  for (const auto j : make_range(dim))
                  {

                    const auto S_ij = S(i,j);
                    const auto S_inv_ji = S_inv(j,i);

                    // Apply the stuff only depending on i and j before entering a, b loops
                    // beta
                    const std::vector<Real> d2beta_dS2_coefs_ij_applied{
                      d2beta_dS2_coefs_times_distortion_weight[1] * S_ij,
                      d2beta_dS2_coefs_times_distortion_weight[2] * S_ij,
                      d2beta_dS2_coefs_times_distortion_weight[3] * S_inv_ji,
                      d2beta_dS2_coefs_times_distortion_weight[4] * S_inv_ji,
                    };
                    // mu
                    const auto dalpha_dS_coef_ij_applied = dalpha_dS_coef_times_dilation_weight * S_inv_ji;

                    Real d2E_dSdR_l = 0.;
                    Real d2E_dSdR_p = 0.;

                    for (const auto a : make_range(i + 1))
                    {

                      // If this condition is met, both the ijab and abij
                      // contributions to the Jacobian are zero due to the
                      // spasity patterns of dS_dR_l and dS_dR_p and this
                      // iteration may be skipped
                      if (!(a == var_id1 && i == var_id2) &&
                          !(a == var_id2 && i == var_id1))
                        continue;

                      const auto S_inv_ja = S_inv(j,a);

                      const Real d2beta_dS2_coef_ia_applied = d2beta_dS2_coefs_times_distortion_weight[0] * I(i,a);
                      const Real d2beta_dS2_coef_ja_applied = d2beta_dS2_coefs_times_distortion_weight[5] * S_inv_ja;
                      const Real alpha_ja_applied = alpha_times_dilation_weight * S_inv_ja;

                      const auto b_limit = (a == i) ? j + 1 : dim;
                      for (const auto b : make_range(b_limit))
                      {

                        // Combine precomputed coefficients with tensor products
                        // to get d2(beta) / dS2
                        Real d2beta_dS2_times_distortion_weight = (
                          d2beta_dS2_coef_ia_applied     * I(j,b) +
                          d2beta_dS2_coefs_ij_applied[0] * S(a,b)     +
                          d2beta_dS2_coefs_ij_applied[1] * S_inv(b,a) +
                          d2beta_dS2_coefs_ij_applied[2] * S(a,b)     +
                          d2beta_dS2_coefs_ij_applied[3] * S_inv(b,a) +
                          d2beta_dS2_coef_ja_applied * S_inv(b,i)
                        );

                        // Incorporate tensor product portion to get d2(mu) /
                        // dS2
                        const Real d2mu_dS2_times_dilation_weight = dalpha_dS_coef_ij_applied * S_inv(b,a) - alpha_ja_applied * S_inv(b,i);


                        // Chain rule to change d/dS to d/dR
                        const auto d2E_dS2 =
                            d2beta_dS2_times_distortion_weight +
                            d2mu_dS2_times_dilation_weight;

                        // if !(a == var_id1 (next line) && i == var_id2
                        // (outside 'a' loop)), dS_dR_l(p) multiplier is zero
                        d2E_dSdR_l += d2E_dS2 * dS_dR_l(a, b);

                        if (!(i == a && j == b))
                          // if !(a == var_id2 (next line) && i == var_id1
                          // (outside 'a' loop)), dS_dR_p(l) multiplier is zero
                          d2E_dSdR_p += d2E_dS2 * dS_dR_p(a, b);

                      } // for b
                    } // for a
                    d2E_dR2 +=
                        d2E_dSdR_l * dS_dR_p(i, j) + d2E_dSdR_p * dS_dR_l(i, j);
                  }// for j
                }// for i, end tensor contraction

                // Jacobian contribution
                const Real jacobian_contribution = quad_weights[qp] * d2E_dR2;
                K[var_id1][var_id2](l, p) += jacobian_contribution;
                // Jacobian is symmetric, add contribution to p,l entry
                // Don't get the diagonal twice!
                if (p < l)
                  // Note the transposition of var_id1 and var_id2 as these are also jacobian indices
                  K[var_id2][var_id1](p, l) += jacobian_contribution;

              }// for var_id2
            }// for p
          }// for var_id1
        }// for l
      }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

void VariationalSmootherSystem::compute_mesh_quality_info()
{
  // If the reference volume has not yet been computed, compute it.
  if (std::abs(_ref_vol) < TOLERANCE * TOLERANCE)
    prepare_for_smoothing();

  std::unique_ptr<DiffContext> con = this->build_context();
  FEMContext & femcontext = cast_ref<FEMContext &>(*con);
  this->init_context(femcontext);

  const auto & mesh = this->get_mesh();
  const auto dim = mesh.mesh_dimension();
  const Real half_dim = 0.5 * dim;
  const auto distortion_weight = 1. - _dilation_weight;

  // Make pre-requests before reinit() for efficiency in
  // --enable-deprecated builds, and to avoid errors in
  // --disable-deprecated builds.
  const auto & fe_map = femcontext.get_element_fe(0)->get_fe_map();
  const auto & JxW = fe_map.get_JxW();
  fe_map.get_dxyzdxi();
  fe_map.get_dxyzdeta();
  fe_map.get_dxyzdzeta();

  MeshQualityInfo info;

  for (const auto * elem : mesh.active_local_element_ptr_range())
    {
      femcontext.pre_fe_reinit(*this, elem);
      femcontext.elem_fe_reinit();

      // Element-integrated quantities
      Real det_S_int = 0.;
      Real beta_int = 0.;
      Real mu_int = 0.;
      Real combined_int = 0.;

      for (const auto qp : index_range(JxW))
        {
          const auto det_SxW = JxW[qp] * _target_jacobian_dets[elem->type()][qp];
          det_S_int += det_SxW;

          // Grab the physical-to-reference mapping Jacobian matrix (i.e., "S") at this qp
          RealTensor S = get_jacobian_at_qp(fe_map, dim, qp);

          // Apply target element transformation
          S *= _target_jacobians[elem->type()][qp];

          // Determinant of S
          const auto det = S.det();
          const auto det_sq = det * det;

          if (det > info.max_qp_det_S)
            info.max_qp_det_S = det;
          else if (det < info.min_qp_det_S)
            info.min_qp_det_S = det;

          if (det < TOLERANCE * TOLERANCE)
            info.mesh_is_tangled = true;

          // trace of S^T * S
          const auto tr = trace(S.transpose() * S, dim);

          // The chi function allows us to handle degenerate elements
          const auto chi = chi_epsilon(det, _epsilon_squared_assembly);

          // distortion
          const Real beta = std::pow(tr / dim, half_dim) / chi;
          beta_int += beta * det_SxW;

          // dilation
          const Real mu = 0.5 * (_ref_vol + det_sq / _ref_vol) / chi;
          mu_int += mu * det_SxW;

          // combined
          const Real E = distortion_weight * beta + _dilation_weight * mu;
          combined_int += E * det_SxW;
        }

      info.total_det_S += det_S_int;
      if (det_S_int > info.max_elem_det_S.second)
        info.max_elem_det_S = std::make_pair(elem->id(), det_S_int);
      else if (det_S_int < info.min_elem_det_S.second)
        info.min_elem_det_S = std::make_pair(elem->id(), det_S_int);

      info.total_distortion += beta_int;
      if (beta_int > info.max_elem_distortion.second)
        info.max_elem_distortion = std::make_pair(elem->id(), beta_int);
      else if (beta_int < info.min_elem_distortion.second)
        info.min_elem_distortion = std::make_pair(elem->id(), beta_int);

      info.total_dilation += mu_int;
      if (mu_int > info.max_elem_dilation.second)
        info.max_elem_dilation = std::make_pair(elem->id(), mu_int);
      else if (mu_int < info.min_elem_dilation.second)
        info.min_elem_dilation = std::make_pair(elem->id(), mu_int);

      info.total_combined += combined_int;
      if (combined_int > info.max_elem_combined.second)
        info.max_elem_combined = std::make_pair(elem->id(), combined_int);
      else if (combined_int < info.min_elem_combined.second)
        info.min_elem_combined = std::make_pair(elem->id(), combined_int);

    } // for elem

  // Get contributions from elements on other processors
  mesh.comm().max(info.max_elem_det_S);
  mesh.comm().min(info.min_elem_det_S);
  mesh.comm().max(info.max_qp_det_S);
  mesh.comm().min(info.min_qp_det_S);
  mesh.comm().sum(info.total_det_S);

  mesh.comm().max(info.max_elem_distortion);
  mesh.comm().min(info.min_elem_distortion);
  mesh.comm().sum(info.total_distortion);

  mesh.comm().max(info.max_elem_dilation);
  mesh.comm().min(info.min_elem_dilation);
  mesh.comm().sum(info.total_dilation);

  mesh.comm().max(info.max_elem_combined);
  mesh.comm().min(info.min_elem_combined);
  mesh.comm().sum(info.total_combined);

  _mesh_info = info;
}

std::pair<std::unique_ptr<Elem>, std::vector<std::unique_ptr<Node>>>
VariationalSmootherSystem::get_target_elem(const ElemType & type)
{
  // Build target element
  auto target_elem = Elem::build(type);

  // Volume of reference element
  const auto ref_vol = target_elem->reference_elem()->volume();

  // Update the nodes of the target element, depending on type
  const Real sqrt_3 = std::sqrt(Real(3));
  std::vector<std::unique_ptr<Node>> owned_nodes;

  const auto type_str = Utility::enum_to_string(type);

  // Elems deriving from Tri
  if (type_str.compare(0, 3, "TRI") == 0)
    {

      // The target element will be an equilateral triangle with area equal to
      // the area of the reference element.

      // Equilateral triangle side length preserving area of the reference element
      const auto side_length = std::sqrt(4. / sqrt_3 * ref_vol);

      // Define the nodal locations of the vertices
      const auto & s = side_length;
      //                                         x        y                  node_id
      owned_nodes.emplace_back(Node::build(Point(0.,      0.),               0));
      owned_nodes.emplace_back(Node::build(Point(s,       0.),               1));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * sqrt_3 * s), 2));

      switch (type)
        {
            case TRI3: {
              // Nothing to do here, vertices already added above
              break;
            }

            case TRI6: {
              // Define the midpoint nodes of the equilateral triangle
              //                                         x         y                   node_id
              owned_nodes.emplace_back(Node::build(Point(0.50 * s, 0.00),              3));
              owned_nodes.emplace_back(Node::build(Point(0.75 * s, 0.25 * sqrt_3 * s), 4));
              owned_nodes.emplace_back(Node::build(Point(0.25 * s, 0.25 * sqrt_3 * s), 5));

              break;
            }

          default:
            libmesh_error_msg("Unsupported triangular element: " << type_str);
            break;
        }
    } // if Tri


  // Set the target_elem equal to the reference elem
  else
    for (const auto & node : target_elem->reference_elem()->node_ref_range())
      owned_nodes.emplace_back(Node::build(node, node.id()));

  // Set nodes of target element
  for (const auto & node_ptr : owned_nodes)
    target_elem->set_node(node_ptr->id(), node_ptr.get());

  libmesh_assert(relative_fuzzy_equals(target_elem->volume(), ref_vol));

  return std::make_pair(std::move(target_elem), std::move(owned_nodes));
}

void VariationalSmootherSystem::get_target_to_reference_jacobian(
    const Elem * const target_elem,
    const FEMContext & femcontext,
    std::vector<RealTensor> & jacobians,
    std::vector<Real> & jacobian_dets)
{

  const auto dim = target_elem->dim();

  const auto & qrule_points = femcontext.get_element_qrule().get_points();
  const auto & qrule_weights = femcontext.get_element_qrule().get_weights();
  const auto nq_points = femcontext.get_element_qrule().n_points();

  // If the target element is the reference element, Jacobian matrix is
  // identity, det of inverse is 1. These will only be overwritten if a
  // different target element is explicitly specified.
  jacobians = std::vector<RealTensor>(nq_points, RealTensor(
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.));
  jacobian_dets = std::vector<Real>(nq_points, 1.0);

  // Don't use "if (*target_elem == *(target_elem->reference_elem()))" here, it
  // only compares global node ids, not the node locations themselves.
  bool target_equals_reference = true;
  const auto * ref_elem = target_elem->reference_elem();
  for (const auto local_id : make_range(target_elem->n_nodes()))
    target_equals_reference &= target_elem->node_ref(local_id) == ref_elem->node_ref(local_id);
  if (target_equals_reference)
    return;

  // Create FEMap to compute target_element mapping information
  FEMap fe_map_target;

  // pre-request mapping derivatives
  fe_map_target.get_dxyzdxi();
  fe_map_target.get_dxyzdeta();
  fe_map_target.get_dxyzdzeta();

  // build map
  fe_map_target.init_reference_to_physical_map(dim, qrule_points, target_elem);
  fe_map_target.compute_map(dim, qrule_weights, target_elem, /*d2phi=*/false);

  for (const auto qp : make_range(nq_points))
    {
      // We use Larisa's H notation to denote the reference-to-target jacobian
      RealTensor H = get_jacobian_at_qp(fe_map_target, dim, qp);

      // The target-to-reference jacobian is the inverse of the
      // reference-to-target jacobian
      jacobians[qp] = H.inverse();
      jacobian_dets[qp] = jacobians[qp].det();
    }
}

} // namespace libMesh
