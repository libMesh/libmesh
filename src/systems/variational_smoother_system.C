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
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_ghost_sync.h"

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
      (*node)(d) = (*current_local_solution)(dof_id);
    }
  }

  SyncNodalPositions sync_object(mesh);
  Parallel::sync_dofobject_data_by_id (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync_object);

  FEMSystem::assembly(get_residual, get_jacobian, apply_heterogeneous_constraints, apply_no_constraints);
}

void VariationalSmootherSystem::init_data ()
{
  auto & mesh = this->get_mesh();
  // Add a variable for each dimension of the mesh
  // "r0" for x, "r1" for y, "r2" for z
  for (const auto & d : make_range(mesh.mesh_dimension()))
    this->add_variable ("r" + std::to_string(d), static_cast<Order>(_fe_order),
                        Utility::string_to_enum<FEFamily>(_fe_family));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Set the current_local_solution to the current mesh
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
      auto & fe_map = my_fe->get_fe_map();

      fe_map.get_dxyzdxi();
      fe_map.get_dxyzdeta();
      fe_map.get_dxyzdzeta();
    }

  // Build a corresponding context for the input system if we haven't
  // already
  auto & input_context = input_contexts[&c];
  if (input_system && !input_context)
    {
      input_context = std::make_unique<FEMContext>(*input_system);
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

  const auto * my_fe = c.get_element_fe(0);

  const auto & fe_map = my_fe->get_fe_map();

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
  const unsigned int n_qpoints = c.get_element_qrule().n_points();
  const auto qpoints = c.get_element_qrule().get_points();
  const auto qweights = c.get_element_qrule().get_weights();

  FEMContext & input_c = *libmesh_map_find(input_contexts, &c);
  if (input_system)
    {
      input_c.pre_fe_reinit(*input_system, &elem);
      input_c.elem_fe_reinit();
    }

  for (const auto qp : make_range(n_qpoints))
    {
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
          ).transpose();
          break;

        default:
          libmesh_error_msg("Unsupported dimension.");
      }

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

      // inverse of S
      const RealTensor S_inv = S.inverse();
      // inverse transpose of S
      const RealTensor S_inv_T = S_inv.transpose();

      // Identity matrix
      const RealTensor I(1, 0, 0,
                         0, 1, 0,
                         0, 0, 1);

      // Distortion metric
      const Real chi = 0.5 * (det + std::sqrt(_epsilon_squared + det_sq));
      const Real chi_sq = chi * chi;
      const Real chi_cube = chi_sq * chi;
      const Real sqrt_term = std::sqrt(_epsilon_squared + det_sq);
      // dchi(x) / dx
      const Real chi_prime = 0.5 * (1. + det / sqrt_term);
      // dchi(det(S) / dS
      RealTensor dchi_dS = chi_prime * det * S_inv_T;
      // d2chi(x) / dx2
      const Real chi_2prime = 0.5 * (1. / sqrt_term - det_sq / std::pow(sqrt_term, 3));

      //const Real beta = std::pow(tr_div_dim, 0.5 * dim) / chi;
      const RealTensor dbeta_dS = (std::pow(tr_div_dim, 0.5 * dim - 1) / chi) * S - (std::pow(tr_div_dim, 0.5 * dim) / chi_sq) * dchi_dS;

      // Dilation metric
      //const Real mu = 0.5 * (_ref_vol + det_sq / _ref_vol) / chi;
      RealTensor dmu_dS = ((det_sq /_ref_vol / chi)) * S_inv_T - (0.5 * (_ref_vol + det_sq / _ref_vol) / chi_sq) * dchi_dS;

      // Combined metric
      //const Real E = (1. - _dilation_weight) * beta + _dilation_weight * mu;
      const RealTensor dE_dS = (1. - _dilation_weight) * dbeta_dS + _dilation_weight * dmu_dS;

      // Factor in nonzero contributions of the gradient of (mapping) jacobian
      std::vector<std::vector<std::vector<Real>>> dphi_maps = {dphidxi_map, dphideta_map, dphidzeta_map};
      for (const auto l : elem.node_index_range())
      {
        for (const auto var_id : make_range(dim))
        {
          RealTensor dS_dR = RealTensor(0);
          for (const auto jj : make_range(dim))
            dS_dR(var_id, jj) = dphi_maps[jj][l][qp];

          F[var_id](l) += qweights[qp] * dE_dS.contract(dS_dR);
        }// for var_id
      }// for l


      if (request_jacobian)
      {
        // Compute jacobian of the smoothing system (i.e., the Hessian)

        // Precompute coefficients to be applied to each component of tensor
        // products in the loops below
        const std::vector<Real> d2beta_dS2_coefs = {
          //Part 1
          //
          // multiplies I[i,k] x I[j,l]
          std::pow(tr_div_dim, 0.5 * dim - 1.) / chi,
          // multiplies S[k,l] x S[i,j]
          ((dim - 2.) / dim) * std::pow(tr_div_dim, 0.5 * dim - 2.) / chi,
          // multiplies S_inv[l,k] * S[i,j] (from  dchi_dS x S)
          -(std::pow(tr_div_dim, 0.5 * dim - 1.) / chi_sq) * chi_prime * det,
          //
          //Part 2
          //
          // multiplies S[k,l] x S_inv[j,i] (from S x dchi_dS)
          -(std::pow(tr_div_dim, 0.5 * dim - 1.) / chi_sq) * chi_prime * det,
          // multiplies S_inv[l,k] x S_inv[j,i]
          2. * std::pow(tr_div_dim, 0.5 * dim) * chi_prime * det / chi_cube * chi_prime * det,
          // multiplies S_inv[l,k] x S_inv[j,i] term from d2chi_dS2
          -(std::pow(tr_div_dim, 0.5 * dim) / chi_sq) * chi_2prime * det_sq,
          // multiplies S_inv[l,k] x S_inv[j,i] term
          -(std::pow(tr_div_dim, 0.5 * dim) / chi_sq) * chi_prime * det,
          // multiplies S_inv[j,k] x S_inv[l,i]
          (std::pow(tr_div_dim, 0.5 * dim) / chi_sq) * chi_prime * det,
        };


        // We represent d mu / dS as alpha(S) * S^-T, where alpha is a scalar function
        const Real alpha = (-chi_prime * det_cube + 2. * det_sq * chi - ref_vol_sq * det) / (2. * _ref_vol * chi_sq);
        // d alpha / dS has the form c(S) * S^-T, where c is the scalar coefficient defined below
        const Real dalpha_dS_coef = det * (-4. * _ref_vol * alpha * chi_prime * chi - chi_2prime * det_cube - chi_prime * det_sq + 4 * chi * det - ref_vol_sq) / (2. * _ref_vol * chi_sq);

        for (const auto l: elem.node_index_range())
        {
          for (const auto var_id1 : make_range(dim))
          {
            RealTensor dS_dR_l = RealTensor(0);
            for (const auto ii : make_range(dim))
              dS_dR_l(var_id1, ii) = dphi_maps[ii][l][qp];

            for (const auto p: elem.node_index_range())
            {
              for (const auto var_id2 : make_range(dim))
              {
                RealTensor dS_dR_p = RealTensor(0);
                for (const auto jj : make_range(dim))
                  dS_dR_p(var_id2, jj) = dphi_maps[jj][p][qp];

                Real d2beta_dR2 = 0.;
                Real d2mu_dR2 = 0.;
                for (const auto i : make_range(dim))
                {
                  for (const auto j : make_range(dim))
                  {
                    for (const auto a : make_range(dim))
                    {
                      for (const auto b : make_range(dim))
                      {
                        const std::vector<Real> d2beta_dS2_tensor_contributions =
                        {
                          I(i,a)     * I(j,b),
                          S(a,b)     * S(i,j),
                          S_inv(b,a) * S(i,j),
                          S(a,b)     * S_inv(j,i),
                          S_inv(b,a) * S_inv(j,i),
                          S_inv(b,a) * S_inv(j,i),
                          S_inv(b,a) * S_inv(j,i),
                          S_inv(j,a) * S_inv(b,i),
                        };

                        Real d2beta_dS2 = 0.;
                        for (const auto comp_id : index_range(d2beta_dS2_coefs))
                        {
                          const Real contribution = d2beta_dS2_coefs[comp_id] * d2beta_dS2_tensor_contributions[comp_id];
                          d2beta_dS2 += contribution;
                        }

                        const Real d2mu_dS2 = dalpha_dS_coef * S_inv(b,a) * S_inv(j,i) - alpha * S_inv(b,i) * S_inv(j,a);

                        // Chain rule to change d/dS to d/dR
                        d2beta_dR2 += d2beta_dS2 * dS_dR_l(a, b) * dS_dR_p(i, j);
                        d2mu_dR2 += d2mu_dS2 * dS_dR_l(a, b) * dS_dR_p(i, j);

                      }// for b
                    }// for a
                  }// for j
                }// for i

                const Real contribution = qweights[qp] * ((1. - _dilation_weight) * d2beta_dR2 + _dilation_weight * d2mu_dR2);
                K[var_id1][var_id2](l, p) += contribution;

              }// for var_id2
            }// for p
          }// for var_id1
        }// for l
      }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

} // namespace libMesh
