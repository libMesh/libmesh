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
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

// C++ includes
#include <functional> // std::reference_wrapper

namespace libMesh
{

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
  std::unique_ptr<DiffContext> con = this->build_context();
  FEMContext & femcontext = cast_ref<FEMContext &>(*con);
  this->init_context(femcontext);

  const auto & mesh = this->get_mesh();
  const auto dim = mesh.mesh_dimension();

  Real elem_averaged_det_J_sum = 0.;

  // Make pre-requests before reinit() for efficiency in
  // --enable-deprecated builds, and to avoid errors in
  // --disable-deprecated builds.
  const auto & fe_map = femcontext.get_element_fe(0)->get_fe_map();
  const auto & JxW = fe_map.get_JxW();

  std::map<ElemType, std::vector<Real>> target_elem_inverse_jacobian_dets;

  for (const auto * elem : mesh.active_local_element_ptr_range())
  {
    femcontext.pre_fe_reinit(*this, elem);
    femcontext.elem_fe_reinit();

    // Add target element info, if applicable
    if (_target_inverse_jacobians.find(elem->type()) == _target_inverse_jacobians.end())
    {
      // Create FEMap to compute target_element mapping information
      FEMap fe_map_target;

      // pre-request mapping derivatives
      const auto & dxyzdxi = fe_map_target.get_dxyzdxi();
      const auto & dxyzdeta = fe_map_target.get_dxyzdeta();
      // const auto & dxyzdzeta = fe_map_target.get_dxyzdzeta();

      const auto & qrule_points = femcontext.get_element_qrule().get_points();
      const auto & qrule_weights = femcontext.get_element_qrule().get_weights();
      const auto nq_points = femcontext.get_element_qrule().n_points();

      // If the target element is the reference element, Jacobian matrix is
      // identity, det of inverse is 1. This will only be overwritten if a
      // different target elemen is explicitly specified.
      target_elem_inverse_jacobian_dets[elem->type()] =
          std::vector<Real>(nq_points, 1.0);

      switch (elem->type())
        {
          case TRI3:
            {
              // Build target element: an equilateral triangle
              Tri3 target_elem;

              // equilateral triangle side length that preserves area of reference
              // element
              const Real sqrt_3 = std::sqrt(3.);
              const auto ref_volume = target_elem.reference_elem()->volume();
              const auto s = std::sqrt(4. / sqrt_3 * ref_volume);

              // Set nodes of target element to form an equilateral triangle
              Node node_0 = Node(0., 0.);
              Node node_1 = Node(s, 0.);
              Node node_2 = Node(0.5 * s, 0.5 * s * sqrt_3);
              target_elem.set_node(0) = &node_0;
              target_elem.set_node(1) = &node_1;
              target_elem.set_node(2) = &node_2;

              // build map
              fe_map_target.init_reference_to_physical_map(dim, qrule_points,
                                                           &target_elem);
              fe_map_target.compute_map(dim, qrule_weights, &target_elem,
                                        /*d2phi=*/false);

              // Yes, triangle-to-triangle mappings have constant Jacobians, but we
              // will keep things general for now
              _target_inverse_jacobians[target_elem.type()] =
                  std::vector<RealTensor>(nq_points);
              for (const auto qp : make_range(nq_points))
              {
                const RealTensor H_inv =
                    RealTensor(dxyzdxi[qp](0), dxyzdeta[qp](0), 0, dxyzdxi[qp](1),
                               dxyzdeta[qp](1), 0, 0, 0, 1)
                        .inverse();

                _target_inverse_jacobians[target_elem.type()][qp] = H_inv;
                target_elem_inverse_jacobian_dets[target_elem.type()][qp] =
                    H_inv.det();
              }

              break;
            }

          default:
            break;
        }
    }

    Real elem_integrated_det_J(0.);
    for (const auto qp : index_range(JxW))
      elem_integrated_det_J +=
          JxW[qp] * target_elem_inverse_jacobian_dets[elem->type()][qp];
    const auto ref_elem_vol = elem->reference_elem()->volume();
    elem_averaged_det_J_sum += elem_integrated_det_J / ref_elem_vol;

  } // for elem

  // Get contributions from elements on other processors
  mesh.comm().sum(elem_averaged_det_J_sum);

  _ref_vol = elem_averaged_det_J_sum / mesh.n_active_elem();
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

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  const auto & fe_map = c.get_element_fe(0)->get_fe_map();

  const auto & dxyzdxi = fe_map.get_dxyzdxi();
  const auto & dxyzdeta = fe_map.get_dxyzdeta();
  const auto & dxyzdzeta = fe_map.get_dxyzdzeta();

  const auto & dphidxi_map = fe_map.get_dphidxi_map();
  const auto & dphideta_map = fe_map.get_dphideta_map();
  const auto & dphidzeta_map = fe_map.get_dphidzeta_map();

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
      RealTensor S;
      // RealTensors are always 3x3, so we will fill any dimensions above dim
      // with 1s on the diagonal. This indicates a 1 to 1 relationship between
      // the physical and reference elements in these extra dimensions.
      switch (dim)
      {
        case 1:
          S = RealTensor(
            dxyzdxi[qp](0), 0, 0,
                         0, 1, 0,
                         0, 0, 1
          );
          break;

        case 2:
          S = RealTensor(
            dxyzdxi[qp](0), dxyzdeta[qp](0), 0,
            dxyzdxi[qp](1), dxyzdeta[qp](1), 0,
            0,              0,               1
          );
          break;

        case 3:
          S = RealTensor(
            dxyzdxi[qp],
            dxyzdeta[qp],
            dxyzdzeta[qp]
          ).transpose();  // Note the transposition!
          break;

        default:
          libmesh_error_msg("Unsupported dimension.");
      }

      // Apply any applicable target element transformations
      if (_target_inverse_jacobians.find(elem.type()) !=
          _target_inverse_jacobians.end())
        S = S * _target_inverse_jacobians[elem.type()][qp];

      // Compute quantities needed for the smoothing algorithm

      // determinant
      const Real det = S.det();
      const Real det_sq = det * det;
      const Real det_cube = det_sq * det;

      const Real ref_vol_sq = _ref_vol * _ref_vol;

      // trace of S^T * S
      // DO NOT USE RealTensor.tr for the trace, it will NOT be correct for
      // 1D and 2D meshes because of our hack of putting 1s in the diagonal of
      // S for the extra dimensions
      const RealTensor ST_S = S.transpose() * S;
      Real tr = 0.0;
      for (const auto i : make_range(dim))
        tr += ST_S(i, i);
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
      // When the element is not degenerate (i.e., det(S) > 0), we can set
      // _epsilon_squared to zero to get chi(det(S)) = det(S)
      const Real epsilon_sq = det > 0. ? 0. : _epsilon_squared;
      const Real chi = 0.5 * (det + std::sqrt(epsilon_sq + det_sq));
      const Real chi_sq = chi * chi;
      const Real sqrt_term = std::sqrt(epsilon_sq + det_sq);
      // dchi(x) / dx
      const Real chi_prime = 0.5 * (1. + det / sqrt_term);
      const Real chi_prime_sq = chi_prime * chi_prime;
      // d2chi(x) / dx2
      const Real chi_2prime = 0.5 * (1. / sqrt_term - det_sq / std::pow(sqrt_term, 3));

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

} // namespace libMesh
