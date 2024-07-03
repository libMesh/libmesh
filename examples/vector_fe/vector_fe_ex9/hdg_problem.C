// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "hdg_problem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/static_condensation.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{
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

HDGProblem::HDGProblem(const Real nu_in, const bool cavity_in)
  : nu(nu_in),
    cavity(cavity_in),
    u_true_soln(nu, cavity),
    v_true_soln(nu, cavity),
    p_true_soln(cavity)
{
}

HDGProblem::~HDGProblem() = default;

void
HDGProblem::init()
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
  global_lm_n_dofs = cavity ? 1 : 0;
}

void
HDGProblem::create_identity_residual(const QBase & quadrature,
                                     const std::vector<Real> & JxW_local,
                                     const std::vector<std::vector<Real>> & phi,
                                     const std::vector<Number> & sol,
                                     const std::size_t n_dofs,
                                     DenseVector<Number> & R)
{
  for (const auto qp : make_range(quadrature.n_points()))
    for (const auto i : make_range(n_dofs))
      R(i) -= JxW_local[qp] * phi[i][qp] * sol[qp];
}

void
HDGProblem::create_identity_jacobian(const QBase & quadrature,
                                     const std::vector<Real> & JxW_local,
                                     const std::vector<std::vector<Real>> & phi,
                                     const std::size_t n_dofs,
                                     DenseMatrix<Number> & J)
{
  for (const auto qp : make_range(quadrature.n_points()))
    for (const auto i : make_range(n_dofs))
      for (const auto j : make_range(n_dofs))
        J(i, j) -= JxW_local[qp] * phi[i][qp] * phi[j][qp];
}

void
HDGProblem::compute_stress(const std::vector<Gradient> & vel_gradient,
                           const unsigned int vel_component,
                           std::vector<Gradient> & sigma)
{
  sigma.resize(qrule->n_points());
  for (const auto qp : make_range(qrule->n_points()))
  {
    Gradient qp_p;
    qp_p(vel_component) = p_sol[qp];
    sigma[qp] = nu * vel_gradient[qp] - qp_p;
  }
}

void
HDGProblem::vector_volume_residual(const std::vector<Gradient> & vector_sol,
                                   const std::vector<Number> & scalar_sol,
                                   DenseVector<Number> & R)
{
  for (const auto qp : make_range(qrule->n_points()))
    for (const auto i : make_range(vector_n_dofs))
    {
      // Vector equation dependence on vector dofs
      R(i) += (*JxW)[qp] * ((*vector_phi)[i][qp] * vector_sol[qp]);

      // Vector equation dependence on scalar dofs
      R(i) += (*JxW)[qp] * ((*div_vector_phi)[i][qp] * scalar_sol[qp]);
    }
}

void
HDGProblem::vector_volume_jacobian(DenseMatrix<Number> & Jqq, DenseMatrix<Number> & Jqs)
{
  for (const auto qp : make_range(qrule->n_points()))
    for (const auto i : make_range(vector_n_dofs))
    {
      // Vector equation dependence on vector dofs
      for (const auto j : make_range(vector_n_dofs))
        Jqq(i, j) += (*JxW)[qp] * ((*vector_phi)[i][qp] * (*vector_phi)[j][qp]);

      // Vector equation dependence on scalar dofs
      for (const auto j : make_range(scalar_n_dofs))
        Jqs(i, j) += (*JxW)[qp] * ((*div_vector_phi)[i][qp] * (*scalar_phi)[j][qp]);
    }
}

RealVectorValue
HDGProblem::vel_cross_vel_residual(const std::vector<Number> & u_sol_local,
                                   const std::vector<Number> & v_sol_local,
                                   const unsigned int qp,
                                   const unsigned int vel_component) const
{
  const RealVectorValue U(u_sol_local[qp], v_sol_local[qp]);
  return U * U(vel_component);
}

RealVectorValue
HDGProblem::vel_cross_vel_jacobian(const std::vector<Number> & u_sol_local,
                                   const std::vector<Number> & v_sol_local,
                                   const unsigned int qp,
                                   const unsigned int vel_component,
                                   const unsigned int vel_j_component,
                                   const std::vector<std::vector<Real>> & phi,
                                   const unsigned int j) const
{
  const RealVectorValue U(u_sol_local[qp], v_sol_local[qp]);
  RealVectorValue vector_phi_local;
  vector_phi_local(vel_j_component) = phi[j][qp];
  auto ret = vector_phi_local * U(vel_component);
  if (vel_component == vel_j_component)
    ret += U * phi[j][qp];
  return ret;
}

void
HDGProblem::scalar_volume_residual(const std::vector<Gradient> & vel_gradient,
                                   const unsigned int vel_component,
                                   std::vector<Gradient> & sigma,
                                   DenseVector<Number> & R)
{
  compute_stress(vel_gradient, vel_component, sigma);
  const auto & mms_info = vel_component == 0 ? static_cast<const ExactSoln &>(u_true_soln)
                                             : static_cast<const ExactSoln &>(v_true_soln);
  for (const auto qp : make_range(qrule->n_points()))
  {
    const auto vel_cross_vel = vel_cross_vel_residual(u_sol, v_sol, qp, vel_component);

    // Prepare forcing function
    Real f = 0;
    if (mms)
      f = mms_info.forcing((*q_point)[qp]);

    for (const auto i : make_range(scalar_n_dofs))
    {
      // Scalar equation dependence on vector and pressure dofs
      R(i) += (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * sigma[qp]);

      // Scalar equation dependence on scalar dofs
      R(i) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vel_cross_vel);

      if (mms)
        // Scalar equation RHS
        R(i) -= (*JxW)[qp] * (*scalar_phi)[i][qp] * f;
    }
  }
}

void
HDGProblem::scalar_volume_jacobian(const unsigned int vel_component,
                                   DenseMatrix<Number> & Jsq,
                                   DenseMatrix<Number> & Jsp,
                                   DenseMatrix<Number> & Jsu,
                                   DenseMatrix<Number> & Jsv)
{
  for (const auto qp : make_range(qrule->n_points()))
    for (const auto i : make_range(scalar_n_dofs))
    {
      // Scalar equation dependence on vector dofs
      for (const auto j : make_range(vector_n_dofs))
        Jsq(i, j) += (*JxW)[qp] * nu * ((*grad_scalar_phi)[i][qp] * (*vector_phi)[j][qp]);

      // Scalar equation dependence on pressure dofs
      for (const auto j : make_range(p_n_dofs))
      {
        Gradient p_phi;
        p_phi(vel_component) = (*scalar_phi)[j][qp];
        Jsp(i, j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * p_phi);
      }

      // Scalar equation dependence on scalar dofs
      for (const auto j : make_range(scalar_n_dofs))
      {
        // derivatives wrt 0th component
        {
          const auto vel_cross_vel =
              vel_cross_vel_jacobian(u_sol, v_sol, qp, vel_component, 0, (*scalar_phi), j);
          Jsu(i, j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vel_cross_vel);
        }
        // derivatives wrt 1th component
        {
          const auto vel_cross_vel =
              vel_cross_vel_jacobian(u_sol, v_sol, qp, vel_component, 1, (*scalar_phi), j);
          Jsv(i, j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vel_cross_vel);
        }
      }
    }
}

void
HDGProblem::pressure_volume_residual(DenseVector<Number> & Rp, DenseVector<Number> & Rglm)
{
  for (const auto qp : make_range(qrule->n_points()))
  {
    // Prepare forcing function
    Real f = 0;
    if (mms)
      f = p_true_soln.forcing((*q_point)[qp]);

    const Gradient vel(u_sol[qp], v_sol[qp]);
    for (const auto i : make_range(p_n_dofs))
    {
      Rp(i) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * vel);

      if (mms)
        // Pressure equation RHS
        Rp(i) -= (*JxW)[qp] * (*scalar_phi)[i][qp] * f;

      if (cavity)
        Rp(i) -= (*JxW)[qp] * (*scalar_phi)[i][qp] * global_lm_dof_value;
    }

    if (cavity)
      Rglm(0) -= (*JxW)[qp] * p_sol[qp];
  }
}

void
HDGProblem::pressure_volume_jacobian(DenseMatrix<Number> & Jpu,
                                     DenseMatrix<Number> & Jpv,
                                     DenseMatrix<Number> & Jpglm,
                                     DenseMatrix<Number> & Jglmp)
{
  for (const auto qp : make_range(qrule->n_points()))
  {
    for (const auto i : make_range(p_n_dofs))
    {
      for (const auto j : make_range(scalar_n_dofs))
      {
        {
          const Gradient phi((*scalar_phi)[j][qp], 0);
          Jpu(i, j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * phi);
        }
        {
          const Gradient phi(0, (*scalar_phi)[j][qp]);
          Jpv(i, j) -= (*JxW)[qp] * ((*grad_scalar_phi)[i][qp] * phi);
        }
      }
      if (cavity)
        Jpglm(i, 0) -= (*JxW)[qp] * (*scalar_phi)[i][qp];
    }

    if (cavity)
    {
      libmesh_assert(scalar_n_dofs == p_n_dofs);
      for (const auto j : make_range(p_n_dofs))
        Jglmp(0, j) -= (*JxW)[qp] * (*scalar_phi)[j][qp];
    }
  }
}

void
HDGProblem::pressure_face_residual(DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    const Gradient vel(lm_u_sol[qp], lm_v_sol[qp]);
    const auto vdotn = vel * (*normals)[qp];
    for (const auto i : make_range(p_n_dofs))
      R(i) += (*JxW_face)[qp] * vdotn * (*scalar_phi_face)[i][qp];
  }
}

void
HDGProblem::pressure_face_jacobian(DenseMatrix<Number> & Jplm_u, DenseMatrix<Number> & Jplm_v)
{
  for (const auto qp : make_range(qface->n_points()))
    for (const auto i : make_range(p_n_dofs))
      for (const auto j : make_range(lm_n_dofs))
      {
        {
          const Gradient phi((*lm_phi_face)[j][qp], 0);
          Jplm_u(i, j) += (*JxW_face)[qp] * phi * (*normals)[qp] * (*scalar_phi_face)[i][qp];
        }
        {
          const Gradient phi(0, (*lm_phi_face)[j][qp]);
          Jplm_v(i, j) += (*JxW_face)[qp] * phi * (*normals)[qp] * (*scalar_phi_face)[i][qp];
        }
      }
}

RealVectorValue
HDGProblem::get_dirichlet_velocity(const unsigned int qp) const
{
  if (mms)
    return RealVectorValue(u_true_soln((*qface_point)[qp]), v_true_soln((*qface_point)[qp]));

  if (cavity)
  {
    if (current_bnd == top_bnd)
      return RealVectorValue(1, 0);
    else
      return RealVectorValue(0, 0);
  }
  else // channel flow case
  {
    libmesh_assert(current_bnd != right_bnd);

    if (current_bnd == left_bnd)
      return RealVectorValue(1, 0);

    else
      return RealVectorValue(0, 0);
  }
}

void
HDGProblem::pressure_dirichlet_residual(DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    const auto dirichlet_velocity = get_dirichlet_velocity(qp);
    const auto vdotn = dirichlet_velocity * (*normals)[qp];
    for (const auto i : make_range(p_n_dofs))
      R(i) += (*JxW_face)[qp] * vdotn * (*scalar_phi_face)[i][qp];
  }
}

void
HDGProblem::vector_dirichlet_residual(const unsigned int vel_component, DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    const auto scalar_value = get_dirichlet_velocity(qp)(vel_component);

    // External boundary -> Dirichlet faces -> Vector equation RHS
    for (const auto i : make_range(vector_n_dofs))
      R(i) -= (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * scalar_value;
  }
}

void
HDGProblem::vector_face_residual(const std::vector<Number> & lm_sol, DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
    // Vector equation dependence on LM dofs
    for (const auto i : make_range(vector_n_dofs))
      R(i) -= (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * lm_sol[qp];
}

void
HDGProblem::vector_face_jacobian(DenseMatrix<Number> & Jqlm)
{
  for (const auto qp : make_range(qface->n_points()))
    // Vector equation dependence on LM dofs
    for (const auto i : make_range(vector_n_dofs))
      for (const auto j : make_range(lm_n_dofs))
        Jqlm(i, j) -=
            (*JxW_face)[qp] * ((*vector_phi_face)[i][qp] * (*normals)[qp]) * (*lm_phi_face)[j][qp];
}

void
HDGProblem::scalar_dirichlet_residual(const std::vector<Gradient> & vector_sol,
                                      const std::vector<Number> & scalar_sol,
                                      const unsigned int vel_component,
                                      DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    Gradient qp_p;
    qp_p(vel_component) = p_sol[qp];

    const auto dirichlet_velocity = get_dirichlet_velocity(qp);
    const auto scalar_value = dirichlet_velocity(vel_component);
    ;

    for (const auto i : make_range(scalar_n_dofs))
    {
      // vector
      R(i) -= (*JxW_face)[qp] * nu * (*scalar_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

      // pressure
      R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

      // scalar from stabilization term
      R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * scalar_sol[qp] * (*normals)[qp] *
              (*normals)[qp];

      // dirichlet lm from stabilization term
      R(i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * scalar_value * (*normals)[qp] *
              (*normals)[qp];

      // dirichlet lm from advection term
      R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (dirichlet_velocity * (*normals)[qp]) *
              scalar_value;
    }
  }
}

void
HDGProblem::scalar_dirichlet_jacobian(const unsigned int vel_component,
                                      DenseMatrix<Number> & Jsq,
                                      DenseMatrix<Number> & Jsp,
                                      DenseMatrix<Number> & Jss)
{
  for (const auto qp : make_range(qface->n_points()))
    for (const auto i : make_range(scalar_n_dofs))
    {
      for (const auto j : make_range(vector_n_dofs))
        Jsq(i, j) -= (*JxW_face)[qp] * nu * (*scalar_phi_face)[i][qp] *
                     ((*vector_phi_face)[j][qp] * (*normals)[qp]);

      for (const auto j : make_range(p_n_dofs))
      {
        Gradient p_phi;
        p_phi(vel_component) = (*scalar_phi_face)[j][qp];
        // pressure
        Jsp(i, j) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
      }

      for (const auto j : make_range(scalar_n_dofs))
        Jss(i, j) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
                     (*normals)[qp] * (*normals)[qp];
    }
}

void
HDGProblem::scalar_face_residual(const std::vector<Gradient> & vector_sol,
                                 const std::vector<Number> & scalar_sol,
                                 const std::vector<Number> & lm_sol,
                                 const unsigned int vel_component,
                                 DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    Gradient qp_p;
    qp_p(vel_component) = p_sol[qp];
    const auto vel_cross_vel = vel_cross_vel_residual(lm_u_sol, lm_v_sol, qp, vel_component);

    for (const auto i : make_range(scalar_n_dofs))
    {
      if (neigh)
      {
        // vector
        R(i) -=
            (*JxW_face)[qp] * nu * (*scalar_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

        // pressure
        R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

        // scalar from stabilization term
        R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * scalar_sol[qp] *
                (*normals)[qp] * (*normals)[qp];

        // lm from stabilization term
        R(i) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * lm_sol[qp] * (*normals)[qp] *
                (*normals)[qp];
      }
      else
        libmesh_assert_msg(current_bnd == right_bnd,
                           "This method should only be called with an outflow boundary. "
                           "Dirichlet boundaries should call a different routine");

      // lm from convection term
      R(i) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
    }
  }
}

void
HDGProblem::scalar_face_jacobian(const unsigned int vel_component,
                                 DenseMatrix<Number> & Jsq,
                                 DenseMatrix<Number> & Jsp,
                                 DenseMatrix<Number> & Jss,
                                 DenseMatrix<Number> & Jslm,
                                 DenseMatrix<Number> & Js_lmu,
                                 DenseMatrix<Number> & Js_lmv)
{
  for (const auto qp : make_range(qface->n_points()))
    for (const auto i : make_range(scalar_n_dofs))
    {
      if (neigh)
      {
        for (const auto j : make_range(vector_n_dofs))
          Jsq(i, j) -= (*JxW_face)[qp] * nu * (*scalar_phi_face)[i][qp] *
                       ((*vector_phi_face)[j][qp] * (*normals)[qp]);

        for (const auto j : make_range(p_n_dofs))
        {
          Gradient p_phi;
          p_phi(vel_component) = (*scalar_phi_face)[j][qp];
          // pressure
          Jsp(i, j) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
        }

        for (const auto j : make_range(scalar_n_dofs))
          Jss(i, j) += (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau *
                       (*scalar_phi_face)[j][qp] * (*normals)[qp] * (*normals)[qp];
      }
      else
        libmesh_assert_msg(current_bnd == right_bnd,
                           "This method should only be called with an outflow boundary. "
                           "Dirichlet boundaries should call a different routine");

      for (const auto j : make_range(lm_n_dofs))
      {
        if (neigh)
          // from stabilization term
          Jslm(i, j) -= (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * tau * (*lm_phi_face)[j][qp] *
                        (*normals)[qp] * (*normals)[qp];

        //
        // from convection term
        //

        // derivatives wrt 0th component
        {
          const auto vel_cross_vel =
              vel_cross_vel_jacobian(lm_u_sol, lm_v_sol, qp, vel_component, 0, (*lm_phi_face), j);
          Js_lmu(i, j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
        }
        // derivatives wrt 1th component
        {
          const auto vel_cross_vel =
              vel_cross_vel_jacobian(lm_u_sol, lm_v_sol, qp, vel_component, 1, (*lm_phi_face), j);
          Js_lmv(i, j) +=
              (*JxW_face)[qp] * (*scalar_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
        }
      }
    }
}

void
HDGProblem::lm_face_residual(const std::vector<Gradient> & vector_sol,
                             const std::vector<Number> & scalar_sol,
                             const std::vector<Number> & lm_sol,
                             const unsigned int vel_component,
                             DenseVector<Number> & R)
{
  for (const auto qp : make_range(qface->n_points()))
  {
    Gradient qp_p;
    qp_p(vel_component) = p_sol[qp];
    const auto vel_cross_vel = vel_cross_vel_residual(lm_u_sol, lm_v_sol, qp, vel_component);

    for (const auto i : make_range(lm_n_dofs))
    {
      // vector
      R(i) -= (*JxW_face)[qp] * nu * (*lm_phi_face)[i][qp] * (vector_sol[qp] * (*normals)[qp]);

      // pressure
      R(i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (qp_p * (*normals)[qp]);

      // scalar from stabilization term
      R(i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * scalar_sol[qp] * (*normals)[qp] *
              (*normals)[qp];

      // lm from stabilization term
      R(i) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * lm_sol[qp] * (*normals)[qp] *
              (*normals)[qp];

      // If we are an internal face we add the convective term. On the outflow boundary we do not
      // zero out the convection term, e.g. we are going to set q + p + tau * (u - u_hat) to zero
      if (neigh)
        // lm from convection term
        R(i) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
      else
        libmesh_assert_msg(current_bnd == right_bnd,
                           "This method should only be called with an outflow boundary. "
                           "Dirichlet boundaries should call a different routine");
    }
  }
}

void
HDGProblem::lm_face_jacobian(const unsigned int vel_component,
                             DenseMatrix<Number> & Jlmq,
                             DenseMatrix<Number> & Jlmp,
                             DenseMatrix<Number> & Jlms,
                             DenseMatrix<Number> & Jlmlm,
                             DenseMatrix<Number> & Jlm_lmu,
                             DenseMatrix<Number> & Jlm_lmv)
{
  for (const auto qp : make_range(qface->n_points()))
    for (const auto i : make_range(lm_n_dofs))
    {
      for (const auto j : make_range(vector_n_dofs))
        Jlmq(i, j) -= (*JxW_face)[qp] * nu * (*lm_phi_face)[i][qp] *
                      ((*vector_phi_face)[j][qp] * (*normals)[qp]);

      for (const auto j : make_range(p_n_dofs))
      {
        Gradient p_phi;
        p_phi(vel_component) = (*scalar_phi_face)[j][qp];
        Jlmp(i, j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * (p_phi * (*normals)[qp]);
      }

      for (const auto j : make_range(scalar_n_dofs))
        Jlms(i, j) += (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * (*scalar_phi_face)[j][qp] *
                      (*normals)[qp] * (*normals)[qp];

      for (const auto j : make_range(lm_n_dofs))
      {
        // from stabilization term
        Jlmlm(i, j) -= (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * tau * (*lm_phi_face)[j][qp] *
                       (*normals)[qp] * (*normals)[qp];
        if (neigh)
        {
          // derivatives wrt 0th component
          {
            const auto vel_cross_vel =
                vel_cross_vel_jacobian(lm_u_sol, lm_v_sol, qp, vel_component, 0, (*lm_phi_face), j);
            Jlm_lmu(i, j) +=
                (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
          }
          // derivatives wrt 1th component
          {
            const auto vel_cross_vel =
                vel_cross_vel_jacobian(lm_u_sol, lm_v_sol, qp, vel_component, 1, (*lm_phi_face), j);
            Jlm_lmv(i, j) +=
                (*JxW_face)[qp] * (*lm_phi_face)[i][qp] * vel_cross_vel * (*normals)[qp];
          }
        }
        else
          libmesh_assert_msg(current_bnd == right_bnd,
                             "This method should only be called with an outflow boundary. "
                             "Dirichlet boundaries should call a different routine");
      }
    }
}

void
HDGProblem::residual(const NumericVector<Number> & X,
                     NumericVector<Number> & R,
                     NonlinearImplicitSystem & S)
{
  R.zero();

  const auto u_num = S.variable_number("vel_x");
  const auto v_num = S.variable_number("vel_y");
  const auto qu_num = S.variable_number("qu");
  const auto qv_num = S.variable_number("qv");
  const auto lm_u_num = S.variable_number("lm_u");
  const auto lm_v_num = S.variable_number("lm_v");
  const auto p_num = S.variable_number("pressure");
  const auto global_lm_num = cavity ? S.variable_number("global_lm") : invalid_uint;

  std::vector<boundary_id_type> boundary_ids;
  const auto & boundary_info = mesh->get_boundary_info();
  DenseVector<Number> Rqu, Rqv, Ru, Rv, Rlm_u, Rlm_v, Rp, Rglm;

  for (const auto & elem : mesh->active_local_element_ptr_range())
  {
    // Retrive our dof indices for all fields
    dof_map->dof_indices(elem, qu_dof_indices, qu_num);
    dof_map->dof_indices(elem, u_dof_indices, u_num);
    dof_map->dof_indices(elem, lm_u_dof_indices, lm_u_num);
    dof_map->dof_indices(elem, qv_dof_indices, qv_num);
    dof_map->dof_indices(elem, v_dof_indices, v_num);
    dof_map->dof_indices(elem, lm_v_dof_indices, lm_v_num);
    dof_map->dof_indices(elem, p_dof_indices, p_num);
    if (cavity)
    {
      dof_map->dof_indices(elem, global_lm_dof_indices, global_lm_num);
      libmesh_assert(global_lm_dof_indices.size() == 1);
    }

    vector_n_dofs = qu_dof_indices.size();
    scalar_n_dofs = u_dof_indices.size();
    lm_n_dofs = lm_u_dof_indices.size();
    p_n_dofs = p_dof_indices.size();
    libmesh_assert(p_n_dofs == scalar_n_dofs);
    Rqu.resize(vector_n_dofs);
    Rqv.resize(vector_n_dofs);
    Ru.resize(scalar_n_dofs);
    Rv.resize(scalar_n_dofs);
    Rlm_u.resize(lm_n_dofs);
    Rlm_v.resize(lm_n_dofs);
    Rp.resize(p_n_dofs);
    Rglm.resize(global_lm_dof_indices.size());

    // Reinit our volume FE objects
    vector_fe->reinit(elem);
    scalar_fe->reinit(elem);

    libmesh_assert_equal_to(vector_n_dofs, vector_phi->size());
    libmesh_assert_equal_to(scalar_n_dofs, scalar_phi->size());

    // Get our local element dof values
    X.get(qu_dof_indices, qu_dof_values);
    X.get(u_dof_indices, u_dof_values);
    X.get(lm_u_dof_indices, lm_u_dof_values);
    X.get(qv_dof_indices, qv_dof_values);
    X.get(v_dof_indices, v_dof_values);
    X.get(lm_v_dof_indices, lm_v_dof_values);
    X.get(p_dof_indices, p_dof_values);
    if (cavity)
      global_lm_dof_value = X(global_lm_dof_indices[0]);

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
    vector_volume_residual(qu_sol, u_sol, Rqu);
    scalar_volume_residual(qu_sol, 0, sigma_u, Ru);

    // qv and v
    vector_volume_residual(qv_sol, v_sol, Rqv);
    scalar_volume_residual(qv_sol, 1, sigma_v, Rv);

    // p
    pressure_volume_residual(Rp, Rglm);

    for (auto side : elem->side_index_range())
    {
      neigh = elem->neighbor_ptr(side);

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

      if (neigh == nullptr)
      {
        boundary_info.boundary_ids(elem, side, boundary_ids);
        libmesh_assert(boundary_ids.size() == 1);
        current_bnd = boundary_ids[0];
        if (cavity || (current_bnd != right_bnd))
        {
          // qu, u
          vector_dirichlet_residual(0, Rqu);
          scalar_dirichlet_residual(qu_sol, u_sol, 0, Ru);

          // qv, v
          vector_dirichlet_residual(1, Rqv);
          scalar_dirichlet_residual(qv_sol, v_sol, 1, Rv);

          // p
          pressure_dirichlet_residual(Rp);

          // Set the LMs on these Dirichlet boundary faces to 0
          create_identity_residual(*qface, *JxW_face, *lm_phi_face, lm_u_sol, lm_n_dofs, Rlm_u);
          create_identity_residual(*qface, *JxW_face, *lm_phi_face, lm_v_sol, lm_n_dofs, Rlm_v);

          continue;
        }
      }

      //
      // if we got here, then we are on an internal face or an outlet face
      //

      // qu, u, lm_u
      vector_face_residual(lm_u_sol, Rqu);
      scalar_face_residual(qu_sol, u_sol, lm_u_sol, 0, Ru);
      lm_face_residual(qu_sol, u_sol, lm_u_sol, 0, Rlm_u);

      // qv, v, lm_v
      vector_face_residual(lm_v_sol, Rqv);
      scalar_face_residual(qv_sol, v_sol, lm_v_sol, 1, Rv);
      lm_face_residual(qv_sol, v_sol, lm_v_sol, 1, Rlm_v);

      // p
      pressure_face_residual(Rp);
    }

    R.add_vector(Rqu, qu_dof_indices);
    R.add_vector(Rqv, qv_dof_indices);
    R.add_vector(Ru, u_dof_indices);
    R.add_vector(Rv, v_dof_indices);
    R.add_vector(Rlm_u, lm_u_dof_indices);
    R.add_vector(Rlm_v, lm_v_dof_indices);
    R.add_vector(Rp, p_dof_indices);
    if (cavity)
      R.add_vector(Rglm, global_lm_dof_indices);
  }
}

void
HDGProblem::jacobian(const NumericVector<Number> & X,
                     SparseMatrix<Number> & /*J*/,
                     NonlinearImplicitSystem & S)
{
  const auto u_num = S.variable_number("vel_x");
  const auto v_num = S.variable_number("vel_y");
  const auto qu_num = S.variable_number("qu");
  const auto qv_num = S.variable_number("qv");
  const auto lm_u_num = S.variable_number("lm_u");
  const auto lm_v_num = S.variable_number("lm_v");
  const auto p_num = S.variable_number("pressure");
  const auto global_lm_num = cavity ? S.variable_number("global_lm") : invalid_uint;

  std::vector<boundary_id_type> boundary_ids;
  const auto & boundary_info = mesh->get_boundary_info();

  DenseMatrix<Number> Jqu_qu, Jqv_qv, Jqu_u, Jqv_v, Ju_qu, Jv_qv, Ju_p, Jv_p, Ju_u, Jv_u, Ju_v,
      Jv_v, Jp_u, Jp_v, Jp_glm, Jglm_p, Jp_lmu, Jp_lmv, Jqu_lmu, Jqv_lmv, Ju_lmu, Jv_lmv, Jv_lmu,
      Ju_lmv, Jlmu_qu, Jlmv_qv, Jlmu_p, Jlmv_p, Jlmu_u, Jlmv_v, Jlmu_lmu, Jlmv_lmu, Jlmu_lmv,
      Jlmv_lmv;

  for (const auto & elem : mesh->active_local_element_ptr_range())
  {
    // Retrive our dof indices for all fields
    dof_map->dof_indices(elem, qu_dof_indices, qu_num);
    dof_map->dof_indices(elem, u_dof_indices, u_num);
    dof_map->dof_indices(elem, lm_u_dof_indices, lm_u_num);
    dof_map->dof_indices(elem, qv_dof_indices, qv_num);
    dof_map->dof_indices(elem, v_dof_indices, v_num);
    dof_map->dof_indices(elem, lm_v_dof_indices, lm_v_num);
    dof_map->dof_indices(elem, p_dof_indices, p_num);
    if (cavity)
    {
      dof_map->dof_indices(elem, global_lm_dof_indices, global_lm_num);
      libmesh_assert(global_lm_dof_indices.size() == 1);
    }

    vector_n_dofs = qu_dof_indices.size();
    scalar_n_dofs = u_dof_indices.size();
    lm_n_dofs = lm_u_dof_indices.size();
    p_n_dofs = p_dof_indices.size();
    libmesh_assert(p_n_dofs == scalar_n_dofs);

    Jqu_qu.resize(qu_dof_indices.size(), qu_dof_indices.size());
    Jqv_qv.resize(qv_dof_indices.size(), qv_dof_indices.size());
    Jqu_u.resize(qu_dof_indices.size(), u_dof_indices.size());
    Jqv_v.resize(qv_dof_indices.size(), v_dof_indices.size());
    Ju_qu.resize(u_dof_indices.size(), qu_dof_indices.size());
    Jv_qv.resize(v_dof_indices.size(), qv_dof_indices.size());
    Ju_p.resize(u_dof_indices.size(), p_dof_indices.size());
    Jv_p.resize(v_dof_indices.size(), p_dof_indices.size());
    Ju_u.resize(u_dof_indices.size(), u_dof_indices.size());
    Jv_u.resize(v_dof_indices.size(), u_dof_indices.size());
    Ju_v.resize(u_dof_indices.size(), v_dof_indices.size());
    Jv_v.resize(v_dof_indices.size(), v_dof_indices.size());
    Jp_u.resize(p_dof_indices.size(), u_dof_indices.size());
    Jp_v.resize(p_dof_indices.size(), v_dof_indices.size());
    Jp_glm.resize(p_dof_indices.size(), global_lm_dof_indices.size());
    Jglm_p.resize(global_lm_dof_indices.size(), p_dof_indices.size());
    Jp_lmu.resize(p_dof_indices.size(), lm_u_dof_indices.size());
    Jp_lmv.resize(p_dof_indices.size(), lm_v_dof_indices.size());
    Jqu_lmu.resize(qu_dof_indices.size(), lm_u_dof_indices.size());
    Jqv_lmv.resize(qv_dof_indices.size(), lm_v_dof_indices.size());
    Ju_lmu.resize(u_dof_indices.size(), lm_u_dof_indices.size());
    Jv_lmv.resize(v_dof_indices.size(), lm_v_dof_indices.size());
    Jv_lmu.resize(v_dof_indices.size(), lm_u_dof_indices.size());
    Ju_lmv.resize(u_dof_indices.size(), lm_v_dof_indices.size());
    Jlmu_qu.resize(lm_u_dof_indices.size(), qu_dof_indices.size());
    Jlmv_qv.resize(lm_v_dof_indices.size(), qv_dof_indices.size());
    Jlmu_p.resize(lm_u_dof_indices.size(), p_dof_indices.size());
    Jlmv_p.resize(lm_v_dof_indices.size(), p_dof_indices.size());
    Jlmu_u.resize(lm_u_dof_indices.size(), u_dof_indices.size());
    Jlmv_v.resize(lm_v_dof_indices.size(), v_dof_indices.size());
    Jlmu_lmu.resize(lm_u_dof_indices.size(), lm_u_dof_indices.size());
    Jlmv_lmu.resize(lm_v_dof_indices.size(), lm_u_dof_indices.size());
    Jlmu_lmv.resize(lm_u_dof_indices.size(), lm_v_dof_indices.size());
    Jlmv_lmv.resize(lm_v_dof_indices.size(), lm_v_dof_indices.size());

    // Reinit our volume FE objects
    vector_fe->reinit(elem);
    scalar_fe->reinit(elem);

    libmesh_assert_equal_to(vector_n_dofs, vector_phi->size());
    libmesh_assert_equal_to(scalar_n_dofs, scalar_phi->size());

    // Get our local element dof values
    X.get(qu_dof_indices, qu_dof_values);
    X.get(u_dof_indices, u_dof_values);
    X.get(lm_u_dof_indices, lm_u_dof_values);
    X.get(qv_dof_indices, qv_dof_values);
    X.get(v_dof_indices, v_dof_values);
    X.get(lm_v_dof_indices, lm_v_dof_values);
    X.get(p_dof_indices, p_dof_values);
    if (cavity)
      global_lm_dof_value = X(global_lm_dof_indices[0]);

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
    vector_volume_jacobian(Jqu_qu, Jqu_u);
    scalar_volume_jacobian(0, Ju_qu, Ju_p, Ju_u, Ju_v);

    // qv and v
    vector_volume_jacobian(Jqv_qv, Jqv_v);
    scalar_volume_jacobian(1, Jv_qv, Jv_p, Jv_u, Jv_v);

    // p
    pressure_volume_jacobian(Jp_u, Jp_v, Jp_glm, Jglm_p);

    for (auto side : elem->side_index_range())
    {
      neigh = elem->neighbor_ptr(side);

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

      if (neigh == nullptr)
      {
        boundary_info.boundary_ids(elem, side, boundary_ids);
        libmesh_assert(boundary_ids.size() == 1);
        current_bnd = boundary_ids[0];
        if (cavity || (current_bnd != right_bnd))
        {
          // qu, u
          scalar_dirichlet_jacobian(0, Ju_qu, Ju_p, Ju_u);

          // qv, v
          scalar_dirichlet_jacobian(1, Jv_qv, Jv_p, Jv_v);

          // Set the LMs on these Dirichlet boundary faces to 0
          create_identity_jacobian(*qface, *JxW_face, *lm_phi_face, lm_n_dofs, Jlmu_lmu);
          create_identity_jacobian(*qface, *JxW_face, *lm_phi_face, lm_n_dofs, Jlmv_lmv);

          continue;
        }
      }

      //
      // if we got here, then we are on an internal face or an outlet face
      //

      // qu, u, lm_u
      vector_face_jacobian(Jqu_lmu);
      scalar_face_jacobian(0, Ju_qu, Ju_p, Ju_u, Ju_lmu, Ju_lmu, Ju_lmv);
      lm_face_jacobian(0, Jlmu_qu, Jlmu_p, Jlmu_u, Jlmu_lmu, Jlmu_lmu, Jlmu_lmv);

      // qv, v, lm_v
      vector_face_jacobian(Jqv_lmv);
      scalar_face_jacobian(1, Jv_qv, Jv_p, Jv_v, Jv_lmv, Jv_lmu, Jv_lmv);
      lm_face_jacobian(1, Jlmv_qv, Jlmv_p, Jlmv_v, Jlmv_lmv, Jlmv_lmu, Jlmv_lmv);

      // p
      pressure_face_jacobian(Jp_lmu, Jp_lmv);
    }

    sc->add_matrix(*elem, qu_num, qu_num, Jqu_qu);
    sc->add_matrix(*elem, qv_num, qv_num, Jqv_qv);
    sc->add_matrix(*elem, qu_num, u_num, Jqu_u);
    sc->add_matrix(*elem, qv_num, v_num, Jqv_v);
    sc->add_matrix(*elem, u_num, qu_num, Ju_qu);
    sc->add_matrix(*elem, v_num, qv_num, Jv_qv);
    sc->add_matrix(*elem, u_num, p_num, Ju_p);
    sc->add_matrix(*elem, v_num, p_num, Jv_p);
    sc->add_matrix(*elem, u_num, u_num, Ju_u);
    sc->add_matrix(*elem, v_num, u_num, Jv_u);
    sc->add_matrix(*elem, u_num, v_num, Ju_v);
    sc->add_matrix(*elem, v_num, v_num, Jv_v);
    sc->add_matrix(*elem, p_num, u_num, Jp_u);
    sc->add_matrix(*elem, p_num, v_num, Jp_v);
    if (global_lm_num != invalid_uint)
    {
      sc->add_matrix(*elem, p_num, global_lm_num, Jp_glm);
      sc->add_matrix(*elem, global_lm_num, p_num, Jglm_p);
    }
    sc->add_matrix(*elem, p_num, lm_u_num, Jp_lmu);
    sc->add_matrix(*elem, p_num, lm_v_num, Jp_lmv);
    sc->add_matrix(*elem, qu_num, lm_u_num, Jqu_lmu);
    sc->add_matrix(*elem, qv_num, lm_v_num, Jqv_lmv);
    sc->add_matrix(*elem, u_num, lm_u_num, Ju_lmu);
    sc->add_matrix(*elem, v_num, lm_v_num, Jv_lmv);
    sc->add_matrix(*elem, v_num, lm_u_num, Jv_lmu);
    sc->add_matrix(*elem, u_num, lm_v_num, Ju_lmv);
    sc->add_matrix(*elem, lm_u_num, qu_num, Jlmu_qu);
    sc->add_matrix(*elem, lm_v_num, qv_num, Jlmv_qv);
    sc->add_matrix(*elem, lm_u_num, p_num, Jlmu_p);
    sc->add_matrix(*elem, lm_v_num, p_num, Jlmv_p);
    sc->add_matrix(*elem, lm_u_num, u_num, Jlmu_u);
    sc->add_matrix(*elem, lm_v_num, v_num, Jlmv_v);
    sc->add_matrix(*elem, lm_u_num, lm_u_num, Jlmu_lmu);
    sc->add_matrix(*elem, lm_v_num, lm_u_num, Jlmv_lmu);
    sc->add_matrix(*elem, lm_u_num, lm_v_num, Jlmu_lmv);
    sc->add_matrix(*elem, lm_v_num, lm_v_num, Jlmv_lmv);
  }
}

} // namespace libMesh
