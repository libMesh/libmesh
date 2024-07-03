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

#ifndef HDG_PROBLEM_H
#define HDG_PROBLEM_H

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/fe.h"
#include "exact_soln.h"

namespace libMesh
{
class MeshBase;
class DofMap;
class StaticCondensation;
class QBase;

class HDGProblem : public NonlinearImplicitSystem::ComputeResidual,
                   public NonlinearImplicitSystem::ComputeJacobian
{
public:
  HDGProblem(const Real nu_in, const bool cavity_in);
  ~HDGProblem();

  // Global data structures
  System * system;
  const MeshBase * mesh;
  const DofMap * dof_map;
  StaticCondensation * sc;

  // FE objects
  std::unique_ptr<FEVectorBase> vector_fe;
  std::unique_ptr<FEBase> scalar_fe;
  std::unique_ptr<QBase> qrule;
  std::unique_ptr<FEVectorBase> vector_fe_face;
  std::unique_ptr<FEBase> scalar_fe_face;
  std::unique_ptr<FEBase> lm_fe_face;
  std::unique_ptr<QBase> qface;

  // boundary information
  boundary_id_type left_bnd;
  boundary_id_type top_bnd;
  boundary_id_type right_bnd;
  boundary_id_type bottom_bnd;

  // The kinematic viscosity
  Real nu;

  // Whether we performing a cavity simulation
  bool cavity;

  // The true solutions
  const USoln u_true_soln;
  const VSoln v_true_soln;
  const PSoln p_true_soln;

  // Whether we are performing an MMS study
  bool mms;

  void init();

  virtual void residual(const NumericVector<Number> & X,
                        NumericVector<Number> & R,
                        NonlinearImplicitSystem & S) override;

  virtual void jacobian(const NumericVector<Number> & X,
                        SparseMatrix<Number> & /*J*/,
                        NonlinearImplicitSystem & S) override;

private:
  void create_identity_residual(const QBase & quadrature,
                                const std::vector<Real> & JxW_local,
                                const std::vector<std::vector<Real>> & phi,
                                const std::vector<Number> & sol,
                                const std::size_t n_dofs,
                                DenseVector<Number> & R);

  void create_identity_jacobian(const QBase & quadrature,
                                const std::vector<Real> & JxW_local,
                                const std::vector<std::vector<Real>> & phi,
                                const std::size_t n_dofs,
                                DenseMatrix<Number> & J);

  void compute_stress(const std::vector<Gradient> & vel_gradient,
                      const unsigned int vel_component,
                      std::vector<Gradient> & sigma);

  void vector_volume_residual(const std::vector<Gradient> & vector_sol,
                              const std::vector<Number> & scalar_sol,
                              DenseVector<Number> & R);

  void vector_volume_jacobian(DenseMatrix<Number> & Jqq, DenseMatrix<Number> & Jqs);

  RealVectorValue vel_cross_vel_residual(const std::vector<Number> & u_sol_local,
                                         const std::vector<Number> & v_sol_local,
                                         const unsigned int qp,
                                         const unsigned int vel_component) const;

  RealVectorValue vel_cross_vel_jacobian(const std::vector<Number> & u_sol_local,
                                         const std::vector<Number> & v_sol_local,
                                         const unsigned int qp,
                                         const unsigned int vel_component,
                                         const unsigned int vel_j_component,
                                         const std::vector<std::vector<Real>> & phi,
                                         const unsigned int j) const;

  void scalar_volume_residual(const std::vector<Gradient> & vel_gradient,
                              const unsigned int vel_component,
                              std::vector<Gradient> & sigma,
                              DenseVector<Number> & R);

  void scalar_volume_jacobian(const unsigned int vel_component,
                              DenseMatrix<Number> & Jsq,
                              DenseMatrix<Number> & Jsp,
                              DenseMatrix<Number> & Jsu,
                              DenseMatrix<Number> & Jsv);

  void pressure_volume_residual(DenseVector<Number> & Rp, DenseVector<Number> & Rglm);

  void pressure_volume_jacobian(DenseMatrix<Number> & Jpu,
                                DenseMatrix<Number> & Jpv,
                                DenseMatrix<Number> & Jpglm,
                                DenseMatrix<Number> & Jglmp);

  void pressure_face_residual(DenseVector<Number> & R);

  void pressure_face_jacobian(DenseMatrix<Number> & Jplm_u, DenseMatrix<Number> & Jplm_v);

  RealVectorValue get_dirichlet_velocity(const unsigned int qp) const;

  void pressure_dirichlet_residual(DenseVector<Number> & R);

  void vector_dirichlet_residual(const unsigned int vel_component, DenseVector<Number> & R);

  void vector_face_residual(const std::vector<Number> & lm_sol, DenseVector<Number> & R);

  void vector_face_jacobian(DenseMatrix<Number> & Jqlm);

  void scalar_dirichlet_residual(const std::vector<Gradient> & vector_sol,
                                 const std::vector<Number> & scalar_sol,
                                 const unsigned int vel_component,
                                 DenseVector<Number> & R);

  void scalar_dirichlet_jacobian(const unsigned int vel_component,
                                 DenseMatrix<Number> & Jsq,
                                 DenseMatrix<Number> & Jsp,
                                 DenseMatrix<Number> & Jss);

  void scalar_face_residual(const std::vector<Gradient> & vector_sol,
                            const std::vector<Number> & scalar_sol,
                            const std::vector<Number> & lm_sol,
                            const unsigned int vel_component,
                            DenseVector<Number> & R);

  void scalar_face_jacobian(const unsigned int vel_component,
                            DenseMatrix<Number> & Jsq,
                            DenseMatrix<Number> & Jsp,
                            DenseMatrix<Number> & Jss,
                            DenseMatrix<Number> & Jslm,
                            DenseMatrix<Number> & Js_lmu,
                            DenseMatrix<Number> & Js_lmv);

  void lm_face_residual(const std::vector<Gradient> & vector_sol,
                        const std::vector<Number> & scalar_sol,
                        const std::vector<Number> & lm_sol,
                        const unsigned int vel_component,
                        DenseVector<Number> & R);

  void lm_face_jacobian(const unsigned int vel_component,
                        DenseMatrix<Number> & Jlmq,
                        DenseMatrix<Number> & Jlmp,
                        DenseMatrix<Number> & Jlms,
                        DenseMatrix<Number> & Jlmlm,
                        DenseMatrix<Number> & Jlm_lmu,
                        DenseMatrix<Number> & Jlm_lmv);

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
  std::vector<dof_id_type> global_lm_dof_indices;

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
  Number global_lm_dof_value;

  // Number of dofs on elem
  std::size_t vector_n_dofs;
  std::size_t scalar_n_dofs;
  std::size_t lm_n_dofs;
  std::size_t p_n_dofs;
  std::size_t global_lm_n_dofs;

  // Our stabilization coefficient
  static constexpr Real tau = 1;

  // The current boundary ID
  boundary_id_type current_bnd;

  // The current element neighbor
  const Elem * neigh;

  // Computed via the formula: \bar{q} = \frac{1}{|\partial K|} \int_{\partial K} q
  std::vector<Number> qbar;
};

} // namespace libMesh

#endif // HDG_PROBLEM_H
