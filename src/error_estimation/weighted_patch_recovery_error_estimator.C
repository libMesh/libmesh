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


// libmesh includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/threads.h"
#include "libmesh/weighted_patch_recovery_error_estimator.h"
#include "libmesh/enum_error_estimator_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/enum_to_string.h"

// C++ includes
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>     // for std::sqrt std::pow std::abs

namespace libMesh
{

ErrorEstimatorType WeightedPatchRecoveryErrorEstimator::type() const
{
  return WEIGHTED_PATCH_RECOVERY;
}



void WeightedPatchRecoveryErrorEstimator::estimate_error (const System & system,
                                                          ErrorVector & error_per_cell,
                                                          const NumericVector<Number> * solution_vector,
                                                          bool)
{
  LOG_SCOPE("estimate_error()", "WeightedPatchRecoveryErrorEstimator");

  // The current mesh
  const MeshBase & mesh = system.get_mesh();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.max_elem_id());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Prepare current_local_solution to localize a non-standard
  // solution vector if necessary
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }

  //------------------------------------------------------------
  // Iterate over all the active elements in the mesh
  // that live on this processor.
  Threads::parallel_for (ConstElemRange(mesh.active_local_elements_begin(),
                                        mesh.active_local_elements_end(),
                                        200),
                         EstimateError(system,
                                       *this,
                                       error_per_cell)
                         );

  // Each processor has now computed the error contributions
  // for its local elements, and error_per_cell contains 0 for all the
  // non-local elements.  Summing the vector will provide the true
  // value for each element, local or remote
  this->reduce_error(error_per_cell, system.comm());

  // If we used a non-standard solution before, now is the time to fix
  // the current_local_solution
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }
}



void WeightedPatchRecoveryErrorEstimator::EstimateError::operator()(const ConstElemRange & range) const
{
  // The current mesh
  const MeshBase & mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();

  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  //------------------------------------------------------------
  // Iterate over all the elements in the range.
  for (const auto & elem : range)
    {
      // We'll need an index into the error vector
      const dof_id_type e_id=elem->id();

      // We are going to build a patch containing the current element
      // and its neighbors on the local processor
      Patch patch(mesh.processor_id());

      // If we are reusing patches and the current element
      // already has an estimate associated with it, move on the
      // next element
      if (this->error_estimator.patch_reuse && error_per_cell[e_id] != 0)
        continue;

      // If we are not reusing patches or haven't built one containing this element, we build one

      // Use user specified patch size and growth strategy
      patch.build_around_element (elem, error_estimator.target_patch_size,
                                  error_estimator.patch_growth_strategy);

      // Declare a new_error_per_cell vector to hold error estimates
      // from each element in this patch, or one estimate if we are
      // not reusing patches since we will only be computing error for
      // one cell
      std::vector<Real> new_error_per_cell(1, 0.);
      if (this->error_estimator.patch_reuse)
        new_error_per_cell.resize(patch.size(), 0.);

      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
        {
          const auto norm_type = error_estimator.error_norm.type(var);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#ifdef DEBUG
          bool is_valid_norm_type =
            norm_type == L2 ||
            norm_type == H1_SEMINORM ||
            norm_type == H2_SEMINORM ||
            norm_type == H1_X_SEMINORM ||
            norm_type == H1_Y_SEMINORM ||
            norm_type == H1_Z_SEMINORM ||
            norm_type == L_INF ||
            norm_type == W1_INF_SEMINORM ||
            norm_type == W2_INF_SEMINORM;
          libmesh_assert (is_valid_norm_type);
#endif // DEBUG
#else
          libmesh_assert (norm_type == L2 ||
                          norm_type == L_INF ||
                          norm_type == H1_SEMINORM ||
                          norm_type == H1_X_SEMINORM ||
                          norm_type == H1_Y_SEMINORM ||
                          norm_type == H1_Z_SEMINORM ||
                          norm_type == W1_INF_SEMINORM);
#endif

#ifdef DEBUG
          if (var > 0)
            {
              // We can't mix L_inf and L_2 norms
              bool is_valid_norm_combo =
                ((norm_type == L2 ||
                  norm_type == H1_SEMINORM ||
                  norm_type == H1_X_SEMINORM ||
                  norm_type == H1_Y_SEMINORM ||
                  norm_type == H1_Z_SEMINORM ||
                  norm_type == H2_SEMINORM) &&
                 (error_estimator.error_norm.type(var-1) == L2 ||
                  error_estimator.error_norm.type(var-1) == H1_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_X_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_Y_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_Z_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H2_SEMINORM)) ||
                ((norm_type == L_INF ||
                  norm_type == W1_INF_SEMINORM ||
                  norm_type == W2_INF_SEMINORM) &&
                 (error_estimator.error_norm.type(var-1) == L_INF ||
                  error_estimator.error_norm.type(var-1) == W1_INF_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == W2_INF_SEMINORM));
              libmesh_assert (is_valid_norm_combo);
            }
#endif // DEBUG

          // Possibly skip this variable
          if (error_estimator.error_norm.weight(var) == 0.0) continue;

          // The type of finite element to use for this variable
          const FEType & fe_type = dof_map.variable_type (var);

          const Order element_order = fe_type.order + elem->p_level();

          // Finite element object for use in this patch
          std::unique_ptr<FEBase> fe (FEBase::build (dim, fe_type));

          // Build an appropriate Gaussian quadrature rule
          std::unique_ptr<QBase> qrule =
            fe_type.default_quadrature_rule(dim, error_estimator._extra_order);

          // Tell the finite element about the quadrature rule.
          fe->attach_quadrature_rule (qrule.get());

          // Get Jacobian values, etc..
          const std::vector<Real> & JxW = fe->get_JxW();
          const std::vector<Point> & q_point = fe->get_xyz();

          // Get whatever phi/dphi/d2phi values we need.  Avoid
          // getting them unless the requested norm is actually going
          // to use them.

          const std::vector<std::vector<Real>> * phi = nullptr;
          // If we're using phi to assert the correct dof_indices
          // vector size later, then we'll need to get_phi whether we
          // plan to use it or not.
#ifdef NDEBUG
          if (norm_type == L2 ||
              norm_type == L_INF)
#endif
            phi = &(fe->get_phi());

          const std::vector<std::vector<RealGradient>> * dphi = nullptr;
          if (norm_type == H1_SEMINORM ||
              norm_type == H1_X_SEMINORM ||
              norm_type == H1_Y_SEMINORM ||
              norm_type == H1_Z_SEMINORM ||
              norm_type == W1_INF_SEMINORM)
            dphi = &(fe->get_dphi());

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          const std::vector<std::vector<RealTensor>> * d2phi = nullptr;
          if (norm_type == H2_SEMINORM ||
              norm_type == W2_INF_SEMINORM)
            d2phi = &(fe->get_d2phi());
#endif

          // global DOF indices
          std::vector<dof_id_type> dof_indices;

          // Compute the appropriate size for the patch projection matrices
          // and vectors;
          unsigned int matsize = element_order + 1;
          if (dim > 1)
            {
              matsize *= (element_order + 2);
              matsize /= 2;
            }
          if (dim > 2)
            {
              matsize *= (element_order + 3);
              matsize /= 3;
            }

          DenseMatrix<Number> Kp(matsize,matsize);
          DenseVector<Number> F,    Fx,     Fy,     Fz,     Fxy,     Fxz,     Fyz;
          DenseVector<Number> Pu_h, Pu_x_h, Pu_y_h, Pu_z_h, Pu_xy_h, Pu_xz_h, Pu_yz_h;
          if (norm_type == L2 ||
              norm_type == L_INF)
            {
              F.resize(matsize); Pu_h.resize(matsize);
            }
          else if (norm_type == H1_SEMINORM ||
                   norm_type == W1_INF_SEMINORM ||
                   norm_type == H2_SEMINORM ||
                   norm_type == W2_INF_SEMINORM)
            {
              Fx.resize(matsize); Pu_x_h.resize(matsize); // stores xx in W2 cases
#if LIBMESH_DIM > 1
              Fy.resize(matsize); Pu_y_h.resize(matsize); // stores yy in W2 cases
#endif
#if LIBMESH_DIM > 2
              Fz.resize(matsize); Pu_z_h.resize(matsize); // stores zz in W2 cases
#endif
            }
          else if (norm_type == H1_X_SEMINORM)
            {
              Fx.resize(matsize); Pu_x_h.resize(matsize); // Only need to compute the x gradient for the x component seminorm
            }
          else if (norm_type == H1_Y_SEMINORM)
            {
              libmesh_assert(LIBMESH_DIM > 1);
              Fy.resize(matsize); Pu_y_h.resize(matsize); // Only need to compute the y gradient for the y component seminorm
            }
          else if (norm_type == H1_Z_SEMINORM)
            {
              libmesh_assert(LIBMESH_DIM > 2);
              Fz.resize(matsize); Pu_z_h.resize(matsize); // Only need to compute the z gradient for the z component seminorm
            }

#if LIBMESH_DIM > 1
          if (norm_type == H2_SEMINORM ||
              norm_type == W2_INF_SEMINORM)
            {
              Fxy.resize(matsize); Pu_xy_h.resize(matsize);
#if LIBMESH_DIM > 2
              Fxz.resize(matsize); Pu_xz_h.resize(matsize);
              Fyz.resize(matsize); Pu_yz_h.resize(matsize);
#endif
            }
#endif


          //------------------------------------------------------
          // Loop over each element in the patch and compute their
          // contribution to the patch gradient projection.
          for (const auto & e_p : patch)
            {
              // Reinitialize the finite element data for this element
              fe->reinit (e_p);

              // Get the global DOF indices for the current variable
              // in the current element
              dof_map.dof_indices (e_p, dof_indices, var);
              libmesh_assert (dof_indices.size() == phi->size());

              const unsigned int n_dofs =
                cast_int<unsigned int>(dof_indices.size());
              const unsigned int n_qp   = qrule->n_points();

              // Compute the weighted projection components from this cell.
              // \int_{Omega_e} \psi_i \psi_j = \int_{Omega_e} w * du_h/dx_k \psi_i
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  // Construct the shape function values for the patch projection
                  std::vector<Real> psi(specpoly(dim, element_order, q_point[qp], matsize));

                  const unsigned int psi_size = cast_int<unsigned int>(psi.size());

                  // Patch matrix contribution
                  const unsigned int m = Kp.m(), n = Kp.n();
                  for (unsigned int i=0; i<m; i++)
                    for (unsigned int j=0; j<n; j++)
                      Kp(i,j) += JxW[qp]*psi[i]*psi[j];

                  if (norm_type == L2 ||
                      norm_type == L_INF)
                    {
                      // Compute the solution on the current patch element
                      // the quadrature point
                      Number u_h = libMesh::zero;

                      for (unsigned int i=0; i<n_dofs; i++)
                        u_h += (*phi)[i][qp]*system.current_solution (dof_indices[i]);

                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        F(i) = JxW[qp]*u_h*psi[i];

                    }
                  else if (norm_type == H1_SEMINORM ||
                           norm_type == W1_INF_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (std::size_t i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));



                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        {
                          Fx(i) += JxW[qp]*grad_u_h(0)*psi[i];
#if LIBMESH_DIM > 1
                          Fy(i) += JxW[qp]*grad_u_h(1)*psi[i];
#endif
#if LIBMESH_DIM > 2
                          Fz(i) += JxW[qp]*grad_u_h(2)*psi[i];
#endif
                        }
                    }
                  else if (norm_type == H1_X_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));



                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        {
                          Fx(i) += JxW[qp]*grad_u_h(0)*psi[i];
                        }
                    }
#if LIBMESH_DIM > 1
                  else if (norm_type == H1_Y_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));



                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        {
                          Fy(i) += JxW[qp]*grad_u_h(1)*psi[i];
                        }
                    }
#endif // LIBMESH_DIM > 1
#if LIBMESH_DIM > 2
                  else if (norm_type == H1_Z_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));



                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        {
                          Fz(i) += JxW[qp]*grad_u_h(2)*psi[i];
                        }
                    }
#endif // LIBMESH_DIM > 2
                  else if (norm_type == H2_SEMINORM ||
                           norm_type == W2_INF_SEMINORM)
                    {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                      // Compute the hessian on the current patch element
                      // at the quadrature point
                      Tensor hess_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        hess_u_h.add_scaled ((*d2phi)[i][qp],
                                             system.current_solution(dof_indices[i]));



                      // Patch RHS contributions
                      for (unsigned int i=0; i != psi_size; i++)
                        {
                          Fx(i)  += JxW[qp]*hess_u_h(0,0)*psi[i];
#if LIBMESH_DIM > 1
                          Fy(i)  += JxW[qp]*hess_u_h(1,1)*psi[i];
                          Fxy(i) += JxW[qp]*hess_u_h(0,1)*psi[i];
#endif
#if LIBMESH_DIM > 2
                          Fz(i)  += JxW[qp]*hess_u_h(2,2)*psi[i];
                          Fxz(i) += JxW[qp]*hess_u_h(0,2)*psi[i];
                          Fyz(i) += JxW[qp]*hess_u_h(1,2)*psi[i];
#endif
                        }
#else
                      libmesh_error_msg("ERROR:  --enable-second-derivatives is required \nfor _sobolev_order == 2!");
#endif
                    }
                  else
                    libmesh_error_msg("Unsupported error norm type == " << Utility::enum_to_string(norm_type));
                } // end quadrature loop
            } // end patch loop



          //--------------------------------------------------
          // Now we have fully assembled the projection system
          // for this patch.  Project the gradient components.
          // MAY NEED TO USE PARTIAL PIVOTING!
          if (norm_type == L2 ||
              norm_type == L_INF)
            {
              Kp.lu_solve(F, Pu_h);
            }
          else if (norm_type == H1_SEMINORM ||
                   norm_type == W1_INF_SEMINORM ||
                   norm_type == H2_SEMINORM ||
                   norm_type == W2_INF_SEMINORM)
            {
              Kp.lu_solve (Fx, Pu_x_h);
#if LIBMESH_DIM > 1
              Kp.lu_solve (Fy, Pu_y_h);
#endif
#if LIBMESH_DIM > 2
              Kp.lu_solve (Fz, Pu_z_h);
#endif
            }
          else if (norm_type == H1_X_SEMINORM)
            {
              Kp.lu_solve (Fx, Pu_x_h);
            }
          else if (norm_type == H1_Y_SEMINORM)
            {
              Kp.lu_solve (Fy, Pu_y_h);
            }
          else if (norm_type == H1_Z_SEMINORM)
            {
              Kp.lu_solve (Fz, Pu_z_h);
            }

#if LIBMESH_DIM > 1
          if (norm_type == H2_SEMINORM ||
              norm_type == W2_INF_SEMINORM)
            {
              Kp.lu_solve(Fxy, Pu_xy_h);
#if LIBMESH_DIM > 2
              Kp.lu_solve(Fxz, Pu_xz_h);
              Kp.lu_solve(Fyz, Pu_yz_h);
#endif
            }
#endif

          // If we are reusing patches, reuse the current patch to loop
          // over all elements in the current patch, otherwise build a new
          // patch containing just the current element and loop over it
          // Note that C++ will not allow patch_re_end to be a const here
          Patch::const_iterator patch_re_it;
          Patch::const_iterator patch_re_end;

          // Declare a new patch
          Patch patch_re(mesh.processor_id());

          if (this->error_estimator.patch_reuse)
            {
              // Just get the iterators from the current patch
              patch_re_it  = patch.begin();
              patch_re_end = patch.end();
            }
          else
            {
              // Use a target patch size of just 0, this will contain
              // just the current element
              patch_re.build_around_element (elem, 0,
                                             error_estimator.patch_growth_strategy);

              // Get the iterators from this newly constructed patch
              patch_re_it = patch_re.begin();
              patch_re_end = patch_re.end();
            }

          // If we are reusing patches, loop over all the elements
          // in the current patch and develop an estimate
          // for all the elements by computing  || w * (P u_h - u_h)|| or ||w *(P grad_u_h - grad_u_h)||
          // or ||w * (P hess_u_h - hess_u_h)|| according to the requested
          // seminorm, otherwise just compute it for the current element

          // Get an FEMContext for this system, this will help us in
          // obtaining the weights from the user code.
          // We don't use full elem_jacobian or subjacobians here.
          FEMContext femcontext(system, nullptr,
                                /* allocate_local_matrices = */ false);
          error_estimator.weight_functions[var]->init_context(femcontext);

          // Loop over every element in the patch
          for (unsigned int e = 0 ; patch_re_it != patch_re_end; ++patch_re_it, ++e)
            {
              // Build the Finite Element for the current element

              // The pth element in the patch
              const Elem * e_p = *patch_re_it;

              // We'll need an index into the error vector for this element
              const dof_id_type e_p_id = e_p->id();

              // Initialize the FEMContext
              femcontext.pre_fe_reinit(system, e_p);

              // We will update the new_error_per_cell vector with element_error if the
              // error_per_cell[e_p_id] entry is non-zero, otherwise update it
              // with 0. i.e. leave it unchanged

              // No need to compute the estimate if we are reusing patches and already have one
              if (this->error_estimator.patch_reuse && error_per_cell[e_p_id] != 0.)
                continue;

              // Reinitialize the finite element data for this element
              fe->reinit (e_p);

              // Get the global DOF indices for the current variable
              // in the current element
              dof_map.dof_indices (e_p, dof_indices, var);
              libmesh_assert (dof_indices.size() == phi->size());

              // The number of dofs for this variable on this element
              const unsigned int n_dofs =
                cast_int<unsigned int>(dof_indices.size());

              // Variable to hold the error on the current element
              Real element_error = 0;

              const Order qorder = fe_type.order + e_p->p_level();

              // A quadrature rule for this element
              QGrid samprule (dim, qorder);

              if (norm_type == W1_INF_SEMINORM ||
                  norm_type == W2_INF_SEMINORM)
                fe->attach_quadrature_rule (&samprule);

              // The number of points we will sample over
              const unsigned int n_sp =
                cast_int<unsigned int>(JxW.size());

              // Loop over every sample point for the current element
              for (unsigned int sp=0; sp<n_sp; sp++)
                {
                  // Compute the solution at the current sample point

                  std::vector<Number> temperr(6,0.0); // x,y,z or xx,yy,zz,xy,xz,yz

                  if (norm_type == L2 ||
                      norm_type == L_INF)
                    {
                      // Compute the value at the current sample point
                      Number u_h = libMesh::zero;

                      for (unsigned int i=0; i<n_dofs; i++)
                        u_h += (*phi)[i][sp]*system.current_solution (dof_indices[i]);

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));
                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[0] += psi[i]*Pu_h(i);
                        }

                      temperr[0] -= u_h;
                    }
                  else if (norm_type == H1_SEMINORM ||
                           norm_type == W1_INF_SEMINORM)
                    {
                      // Compute the gradient at the current sample point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][sp],
                                             system.current_solution(dof_indices[i]));

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));

                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[0] += psi[i]*Pu_x_h(i);
#if LIBMESH_DIM > 1
                          temperr[1] += psi[i]*Pu_y_h(i);
#endif
#if LIBMESH_DIM > 2
                          temperr[2] += psi[i]*Pu_z_h(i);
#endif
                        }
                      temperr[0] -= grad_u_h(0);
#if LIBMESH_DIM > 1
                      temperr[1] -= grad_u_h(1);
#endif
#if LIBMESH_DIM > 2
                      temperr[2] -= grad_u_h(2);
#endif
                    }
                  else if (norm_type == H1_X_SEMINORM)
                    {
                      // Compute the gradient at the current sample point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][sp],
                                             system.current_solution(dof_indices[i]));

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));
                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[0] += psi[i]*Pu_x_h(i);
                        }

                      temperr[0] -= grad_u_h(0);
                    }
#if LIBMESH_DIM > 1
                  else if (norm_type == H1_Y_SEMINORM)
                    {
                      // Compute the gradient at the current sample point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][sp],
                                             system.current_solution(dof_indices[i]));

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));
                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[1] += psi[i]*Pu_y_h(i);
                        }

                      temperr[1] -= grad_u_h(1);
                    }
#endif // LIBMESH_DIM > 1
#if LIBMESH_DIM > 2
                  else if (norm_type == H1_Z_SEMINORM)
                    {
                      // Compute the gradient at the current sample point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][sp],
                                             system.current_solution(dof_indices[i]));

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));
                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[2] += psi[i]*Pu_z_h(i);
                        }

                      temperr[2] -= grad_u_h(2);
                    }
#endif // LIBMESH_DIM > 2
                  else if (norm_type == H2_SEMINORM ||
                           norm_type == W2_INF_SEMINORM)
                    {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                      // Compute the Hessian at the current sample point
                      Tensor hess_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        hess_u_h.add_scaled ((*d2phi)[i][sp],
                                             system.current_solution(dof_indices[i]));

                      // Compute the phi values at the current sample point
                      std::vector<Real> psi(specpoly(dim, element_order, q_point[sp], matsize));
                      for (unsigned int i=0; i<matsize; i++)
                        {
                          temperr[0] += psi[i]*Pu_x_h(i);
#if LIBMESH_DIM > 1
                          temperr[1] += psi[i]*Pu_y_h(i);
                          temperr[3] += psi[i]*Pu_xy_h(i);
#endif
#if LIBMESH_DIM > 2
                          temperr[2] += psi[i]*Pu_z_h(i);
                          temperr[4] += psi[i]*Pu_xz_h(i);
                          temperr[5] += psi[i]*Pu_yz_h(i);
#endif
                        }

                      temperr[0] -= hess_u_h(0,0);
#if LIBMESH_DIM > 1
                      temperr[1] -= hess_u_h(1,1);
                      temperr[3] -= hess_u_h(0,1);
#endif
#if LIBMESH_DIM > 2
                      temperr[2] -= hess_u_h(2,2);
                      temperr[4] -= hess_u_h(0,2);
                      temperr[5] -= hess_u_h(1,2);
#endif
#else
                      libmesh_error_msg("ERROR:  --enable-second-derivatives is required \nfor _sobolev_order == 2!");
#endif
                    }

                  // Get the weight from the user code
                  Number weight = (*error_estimator.weight_functions[var])(femcontext, q_point[sp], system.time);

                  // Add up relevant terms.  We can easily optimize the
                  // LIBMESH_DIM < 3 cases a little bit with the exception
                  // of the W2 cases

                  if (norm_type == L_INF)
                    element_error = std::max(element_error, std::abs(weight*temperr[0]));
                  else if (norm_type == W1_INF_SEMINORM)
                    for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                      element_error = std::max(element_error, std::abs(weight*temperr[i]));
                  else if (norm_type == W2_INF_SEMINORM)
                    for (unsigned int i=0; i != 6; ++i)
                      element_error = std::max(element_error, std::abs(weight*temperr[i]));
                  else if (norm_type == L2)
                    element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[0]);
                  else if (norm_type == H1_SEMINORM)
                    for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                      element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[i]);
                  else if (norm_type == H1_X_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[0]);
                  else if (norm_type == H1_Y_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[1]);
                  else if (norm_type == H1_Z_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[2]);
                  else if (norm_type == H2_SEMINORM)
                    {
                      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                        element_error += JxW[sp]*TensorTools::norm_sq(weight*temperr[i]);
                      // Off diagonal terms enter into the Hessian norm twice
                      for (unsigned int i=3; i != 6; ++i)
                        element_error += JxW[sp]*2*TensorTools::norm_sq(weight*temperr[i]);
                    }

                } // End loop over sample points

              if (norm_type == L_INF ||
                  norm_type == W1_INF_SEMINORM ||
                  norm_type == W2_INF_SEMINORM)
                new_error_per_cell[e] += error_estimator.error_norm.weight(var) * element_error;
              else if (norm_type == L2 ||
                       norm_type == H1_SEMINORM ||
                       norm_type == H1_X_SEMINORM ||
                       norm_type == H1_Y_SEMINORM ||
                       norm_type == H1_Z_SEMINORM ||
                       norm_type == H2_SEMINORM)
                new_error_per_cell[e] += error_estimator.error_norm.weight_sq(var) * element_error;
              else
                libmesh_error_msg("Unsupported error norm type == " << Utility::enum_to_string(norm_type));
            }  // End (re) loop over patch elements

        } // end variables loop

      // Now that we have the contributions from each variable,
      // we have take square roots of the entries we
      // added to error_per_cell to get an error norm
      // If we are reusing patches, once again reuse the current patch to loop
      // over all elements in the current patch, otherwise build a new
      // patch containing just the current element and loop over it
      Patch::const_iterator patch_re_it;
      Patch::const_iterator patch_re_end;

      // Build a new patch if necessary
      Patch current_elem_patch(mesh.processor_id());

      if (this->error_estimator.patch_reuse)
        {
          // Just get the iterators from the current patch
          patch_re_it  = patch.begin();
          patch_re_end = patch.end();
        }
      else
        {
          // Use a target patch size of just 0, this will contain
          // just the current element.
          current_elem_patch.build_around_element (elem, 0,
                                                   error_estimator.patch_growth_strategy);

          // Get the iterators from this newly constructed patch
          patch_re_it = current_elem_patch.begin();
          patch_re_end = current_elem_patch.end();
        }

      // Loop over every element in the patch we just constructed
      for (unsigned int i = 0 ; patch_re_it != patch_re_end; ++patch_re_it, ++i)
        {
          // The pth element in the patch
          const Elem * e_p = *patch_re_it;

          // We'll need an index into the error vector
          const dof_id_type e_p_id = e_p->id();

          // Update the error_per_cell vector for this element
          if (error_estimator.error_norm.type(0) == L2 ||
              error_estimator.error_norm.type(0) == H1_SEMINORM ||
              error_estimator.error_norm.type(0) == H1_X_SEMINORM ||
              error_estimator.error_norm.type(0) == H1_Y_SEMINORM ||
              error_estimator.error_norm.type(0) == H1_Z_SEMINORM ||
              error_estimator.error_norm.type(0) == H2_SEMINORM)
            {
              Threads::spin_mutex::scoped_lock acquire(Threads::spin_mtx);
              if (!error_per_cell[e_p_id])
                error_per_cell[e_p_id] = static_cast<ErrorVectorReal>
                  (std::sqrt(new_error_per_cell[i]));
            }
          else
            {
              libmesh_assert (error_estimator.error_norm.type(0) == L_INF ||
                              error_estimator.error_norm.type(0) == W1_INF_SEMINORM ||
                              error_estimator.error_norm.type(0) == W2_INF_SEMINORM);
              Threads::spin_mutex::scoped_lock acquire(Threads::spin_mtx);
              if (!error_per_cell[e_p_id])
                error_per_cell[e_p_id] = static_cast<ErrorVectorReal>
                  (new_error_per_cell[i]);
            }

        } // End loop over every element in patch


    } // end element loop

} // End () operator definition

} // namespace libMesh
