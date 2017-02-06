// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>     // for std::sqrt std::pow std::abs


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/error_vector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/patch.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/tensor_value.h"
#include "libmesh/threads.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{
// Setter function for the patch_reuse flag
void PatchRecoveryErrorEstimator::set_patch_reuse(bool patch_reuse_flag)
{
  patch_reuse = patch_reuse_flag;
}

//-----------------------------------------------------------------
// PatchRecoveryErrorEstimator implementations
std::vector<Real> PatchRecoveryErrorEstimator::specpoly(const unsigned int dim,
                                                        const Order order,
                                                        const Point p,
                                                        const unsigned int matsize)
{
  std::vector<Real> psi;
  psi.reserve(matsize);
  unsigned int npows = order+1;
  std::vector<Real> xpow(npows,1.), ypow, zpow;
  {
    Real x = p(0);
    for (unsigned int i=1; i != npows; ++i)
      xpow[i] = xpow[i-1] * x;
  }
  if (dim > 1)
    {
      Real y = p(1);
      ypow.resize(npows,1.);
      for (unsigned int i=1; i != npows; ++i)
        ypow[i] = ypow[i-1] * y;
    }
  if (dim > 2)
    {
      Real z = p(2);
      zpow.resize(npows,1.);
      for (unsigned int i=1; i != npows; ++i)
        zpow[i] = zpow[i-1] * z;
    }

  // builds psi vector of form 1 x y z x^2 xy xz y^2 yz z^2 etc..
  // I haven't added 1D support here
  for (unsigned int poly_deg=0; poly_deg <= static_cast<unsigned int>(order) ; poly_deg++)
    { // loop over all polynomials of total degreee = poly_deg

      switch (dim)
        {
          // 3D spectral polynomial basis functions
        case 3:
          {
            for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
              for (int yexp=poly_deg-xexp; yexp >= 0; yexp--)
                {
                  int zexp = poly_deg - xexp - yexp;
                  psi.push_back(xpow[xexp]*ypow[yexp]*zpow[zexp]);
                }
            break;
          }

          // 2D spectral polynomial basis functions
        case 2:
          {
            for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
              {
                int yexp = poly_deg - xexp;
                psi.push_back(xpow[xexp]*ypow[yexp]);
              }
            break;
          }

          // 1D spectral polynomial basis functions
        case 1:
          {
            int xexp = poly_deg;
            psi.push_back(xpow[xexp]);
            break;
          }

        default:
          libmesh_error_msg("Invalid dimension dim " << dim);
        }
    }

  return psi;
}



void PatchRecoveryErrorEstimator::estimate_error (const System & system,
                                                  ErrorVector & error_per_cell,
                                                  const NumericVector<Number> * solution_vector,
                                                  bool)
{
  LOG_SCOPE("estimate_error()", "PatchRecoveryErrorEstimator");

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



void PatchRecoveryErrorEstimator::EstimateError::operator()(const ConstElemRange & range) const
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
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it!=range.end(); ++elem_it)
    {
      // elem is necessarily an active element on the local processor
      const Elem * elem = *elem_it;

      // We'll need an index into the error vector
      const dof_id_type e_id=elem->id();

      // We are going to build a patch containing the current element
      // and its neighbors on the local processor
      Patch patch(mesh.processor_id());

      // If we are reusing patches and the current element
      // already has an estimate associated with it, move on the
      // next element
      if(this->error_estimator.patch_reuse && error_per_cell[e_id] != 0)
        continue;

      // If we are not reusing patches or havent built one containing this element, we build one

      // Use user specified patch size and growth strategy
      patch.build_around_element (elem, error_estimator.target_patch_size,
                                  error_estimator.patch_growth_strategy);

      // Declare a new_error_per_cell vector to hold error estimates
      // from each element in this patch, or one estimate if we are
      // not reusing patches since we will only be computing error for
      // one cell
      std::vector<Real> new_error_per_cell(1, 0.);
      if(this->error_estimator.patch_reuse)
        new_error_per_cell.resize(patch.size(), 0.);

      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
        {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#ifdef DEBUG
          bool is_valid_norm_type =
            error_estimator.error_norm.type(var) == L2 ||
            error_estimator.error_norm.type(var) == H1_SEMINORM ||
            error_estimator.error_norm.type(var) == H2_SEMINORM ||
            error_estimator.error_norm.type(var) == H1_X_SEMINORM ||
            error_estimator.error_norm.type(var) == H1_Y_SEMINORM ||
            error_estimator.error_norm.type(var) == H1_Z_SEMINORM ||
            error_estimator.error_norm.type(var) == L_INF ||
            error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
            error_estimator.error_norm.type(var) == W2_INF_SEMINORM;
          libmesh_assert (is_valid_norm_type);
#endif // DEBUG
#else
          libmesh_assert (error_estimator.error_norm.type(var) == L2 ||
                          error_estimator.error_norm.type(var) == L_INF ||
                          error_estimator.error_norm.type(var) == H1_SEMINORM ||
                          error_estimator.error_norm.type(var) == H1_X_SEMINORM ||
                          error_estimator.error_norm.type(var) == H1_Y_SEMINORM ||
                          error_estimator.error_norm.type(var) == H1_Z_SEMINORM ||
                          error_estimator.error_norm.type(var) == W1_INF_SEMINORM);
#endif


#ifdef DEBUG
          if (var > 0)
            {
              // We can't mix L_inf and L_2 norms
              bool is_valid_norm_type =
                ((error_estimator.error_norm.type(var) == L2 ||
                  error_estimator.error_norm.type(var) == H1_SEMINORM ||
                  error_estimator.error_norm.type(var) == H1_X_SEMINORM ||
                  error_estimator.error_norm.type(var) == H1_Y_SEMINORM ||
                  error_estimator.error_norm.type(var) == H1_Z_SEMINORM ||
                  error_estimator.error_norm.type(var) == H2_SEMINORM) &&
                 (error_estimator.error_norm.type(var-1) == L2 ||
                  error_estimator.error_norm.type(var-1) == H1_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_X_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_Y_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H1_Z_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == H2_SEMINORM)) ||
                ((error_estimator.error_norm.type(var) == L_INF ||
                  error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                  error_estimator.error_norm.type(var) == W2_INF_SEMINORM) &&
                 (error_estimator.error_norm.type(var-1) == L_INF ||
                  error_estimator.error_norm.type(var-1) == W1_INF_SEMINORM ||
                  error_estimator.error_norm.type(var-1) == W2_INF_SEMINORM));
              libmesh_assert (is_valid_norm_type);
            }
#endif // DEBUG

          // Possibly skip this variable
          if (error_estimator.error_norm.weight(var) == 0.0) continue;

          // The type of finite element to use for this variable
          const FEType & fe_type = dof_map.variable_type (var);

          const Order element_order  = static_cast<Order>
            (fe_type.order + elem->p_level());

          // Finite element object for use in this patch
          UniquePtr<FEBase> fe (FEBase::build (dim, fe_type));

          // Build an appropriate Gaussian quadrature rule
          UniquePtr<QBase> qrule (fe_type.default_quadrature_rule(dim));

          // Tell the finite element about the quadrature rule.
          fe->attach_quadrature_rule (qrule.get());

          // Get Jacobian values, etc..
          const std::vector<Real> & JxW = fe->get_JxW();
          const std::vector<Point> & q_point = fe->get_xyz();

          // Get whatever phi/dphi/d2phi values we need.  Avoid
          // getting them unless the requested norm is actually going
          // to use them.

          const std::vector<std::vector<Real> > * phi = libmesh_nullptr;
          // If we're using phi to assert the correct dof_indices
          // vector size later, then we'll need to get_phi whether we
          // plan to use it or not.
#ifdef NDEBUG
          if (error_estimator.error_norm.type(var) == L2 ||
              error_estimator.error_norm.type(var) == L_INF)
#endif
            phi = &(fe->get_phi());

          const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;
          if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
              error_estimator.error_norm.type(var) == H1_X_SEMINORM ||
              error_estimator.error_norm.type(var) == H1_Y_SEMINORM ||
              error_estimator.error_norm.type(var) == H1_Z_SEMINORM ||
              error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
            dphi = &(fe->get_dphi());

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          const std::vector<std::vector<RealTensor> > * d2phi = libmesh_nullptr;
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            d2phi = &(fe->get_d2phi());
#endif

          // global DOF indices
          std::vector<dof_id_type> dof_indices;

          // Compute the approprite size for the patch projection matrices
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
          if (error_estimator.error_norm.type(var) == L2 ||
              error_estimator.error_norm.type(var) == L_INF)
            {
              F.resize(matsize); Pu_h.resize(matsize);
            }
          else if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
                   error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                   error_estimator.error_norm.type(var) == H2_SEMINORM ||
                   error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
              Fx.resize(matsize); Pu_x_h.resize(matsize); // stores xx in W2 cases
#if LIBMESH_DIM > 1
              Fy.resize(matsize); Pu_y_h.resize(matsize); // stores yy in W2 cases
#endif
#if LIBMESH_DIM > 2
              Fz.resize(matsize); Pu_z_h.resize(matsize); // stores zz in W2 cases
#endif
            }
          else if (error_estimator.error_norm.type(var) == H1_X_SEMINORM)
            {
              Fx.resize(matsize); Pu_x_h.resize(matsize); // Only need to compute the x gradient for the x component seminorm
            }
          else if (error_estimator.error_norm.type(var) == H1_Y_SEMINORM)
            {
              libmesh_assert_greater (LIBMESH_DIM, 1);
              Fy.resize(matsize); Pu_y_h.resize(matsize); // Only need to compute the y gradient for the y component seminorm
            }
          else if (error_estimator.error_norm.type(var) == H1_Z_SEMINORM)
            {
              libmesh_assert_greater (LIBMESH_DIM, 2);
              Fz.resize(matsize); Pu_z_h.resize(matsize); // Only need to compute the z gradient for the z component seminorm
            }

#if LIBMESH_DIM > 1
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
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
          Patch::const_iterator        patch_it  = patch.begin();
          const Patch::const_iterator  patch_end = patch.end();

          for (; patch_it != patch_end; ++patch_it)
            {
              // The pth element in the patch
              const Elem * e_p = *patch_it;

              // Reinitialize the finite element data for this element
              fe->reinit (e_p);

              // Get the global DOF indices for the current variable
              // in the current element
              dof_map.dof_indices (e_p, dof_indices, var);
              libmesh_assert_equal_to (dof_indices.size(), phi->size());

              const unsigned int n_dofs =
                cast_int<unsigned int>(dof_indices.size());
              const unsigned int n_qp   = qrule->n_points();

              // Compute the projection components from this cell.
              // \int_{Omega_e} \psi_i \psi_j = \int_{Omega_e} du_h/dx_k \psi_i
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  // Construct the shape function values for the patch projection
                  std::vector<Real> psi(specpoly(dim, element_order, q_point[qp], matsize));

                  // Patch matrix contribution
                  for (unsigned int i=0; i<Kp.m(); i++)
                    for (unsigned int j=0; j<Kp.n(); j++)
                      Kp(i,j) += JxW[qp]*psi[i]*psi[j];

                  if (error_estimator.error_norm.type(var) == L2 ||
                      error_estimator.error_norm.type(var) == L_INF)
                    {
                      // Compute the solution on the current patch element
                      // the quadrature point
                      Number u_h = libMesh::zero;

                      for (unsigned int i=0; i<n_dofs; i++)
                        u_h += (*phi)[i][qp]*system.current_solution (dof_indices[i]);

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
                        F(i) += JxW[qp]*u_h*psi[i];

                    }
                  else if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
                           error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
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
                  else if (error_estimator.error_norm.type(var) == H1_X_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
                        {
                          Fx(i) += JxW[qp]*grad_u_h(0)*psi[i];
                        }
                    }
                  else if (error_estimator.error_norm.type(var) == H1_Y_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
                        {
                          Fy(i) += JxW[qp]*grad_u_h(1)*psi[i];
                        }
                    }
                  else if (error_estimator.error_norm.type(var) == H1_Z_SEMINORM)
                    {
                      // Compute the gradient on the current patch element
                      // at the quadrature point
                      Gradient grad_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        grad_u_h.add_scaled ((*dphi)[i][qp],
                                             system.current_solution(dof_indices[i]));

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
                        {
                          Fz(i) += JxW[qp]*grad_u_h(2)*psi[i];
                        }
                    }
                  else if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
                           error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
                    {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                      // Compute the hessian on the current patch element
                      // at the quadrature point
                      Tensor hess_u_h;

                      for (unsigned int i=0; i<n_dofs; i++)
                        hess_u_h.add_scaled ((*d2phi)[i][qp],
                                             system.current_solution(dof_indices[i]));

                      // Patch RHS contributions
                      for (std::size_t i=0; i<psi.size(); i++)
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
                    libmesh_error_msg("Unsupported error norm type!");
                } // end quadrature loop
            } // end patch loop



          //--------------------------------------------------
          // Now we have fully assembled the projection system
          // for this patch.  Project the gradient components.
          // MAY NEED TO USE PARTIAL PIVOTING!
          if (error_estimator.error_norm.type(var) == L2 ||
              error_estimator.error_norm.type(var) == L_INF)
            {
              Kp.lu_solve(F, Pu_h);
            }
          else if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
                   error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                   error_estimator.error_norm.type(var) == H2_SEMINORM ||
                   error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
              Kp.lu_solve (Fx, Pu_x_h);
#if LIBMESH_DIM > 1
              Kp.lu_solve (Fy, Pu_y_h);
#endif
#if LIBMESH_DIM > 2
              Kp.lu_solve (Fz, Pu_z_h);
#endif
            }
          else if (error_estimator.error_norm.type(var) == H1_X_SEMINORM)
            {
              Kp.lu_solve (Fx, Pu_x_h);
            }
          else if (error_estimator.error_norm.type(var) == H1_Y_SEMINORM)
            {
              Kp.lu_solve (Fy, Pu_y_h);
            }
          else if (error_estimator.error_norm.type(var) == H1_Z_SEMINORM)
            {
              Kp.lu_solve (Fz, Pu_z_h);
            }

#if LIBMESH_DIM > 1
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
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

          if(this->error_estimator.patch_reuse)
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
          // for all the elements by computing  ||P u_h - u_h|| or ||P grad_u_h - grad_u_h||
          // or ||P hess_u_h - hess_u_h|| according to the requested
          // seminorm, otherwise just compute it for the current element

          // Loop over every element in the patch
          for (unsigned int e = 0 ; patch_re_it != patch_re_end; ++patch_re_it, ++e)
            {
              // Build the Finite Element for the current element

              // The pth element in the patch
              const Elem * e_p = *patch_re_it;

              // We'll need an index into the error vector for this element
              const dof_id_type e_p_id = e_p->id();

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
              libmesh_assert_equal_to (dof_indices.size(), phi->size());

              // The number of dofs for this variable on this element
              const unsigned int n_dofs =
                cast_int<unsigned int>(dof_indices.size());

              // Variable to hold the error on the current element
              Real element_error = 0;

              const Order qorder =
                static_cast<Order>(fe_type.order + e_p->p_level());

              // A quadrature rule for this element
              QGrid samprule (dim, qorder);

              if (error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                  error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
                fe->attach_quadrature_rule (&samprule);

              // The number of points we will sample over
              const unsigned int n_sp =
                cast_int<unsigned int>(JxW.size());

              // Loop over every sample point for the current element
              for (unsigned int sp=0; sp<n_sp; sp++)
                {
                  // Compute the solution at the current sample point

                  std::vector<Number> temperr(6,0.0); // x,y,z or xx,yy,zz,xy,xz,yz

                  if (error_estimator.error_norm.type(var) == L2 ||
                      error_estimator.error_norm.type(var) == L_INF)
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
                  else if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
                           error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
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
                  else if (error_estimator.error_norm.type(var) == H1_X_SEMINORM)
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
                  else if (error_estimator.error_norm.type(var) == H1_Y_SEMINORM)
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
                  else if (error_estimator.error_norm.type(var) == H1_Z_SEMINORM)
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
                  else if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
                           error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
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
                  // Add up relevant terms.  We can easily optimize the
                  // LIBMESH_DIM < 3 cases a little bit with the exception
                  // of the W2 cases

                  if (error_estimator.error_norm.type(var) == L_INF)
                    element_error = std::max(element_error, std::abs(temperr[0]));
                  else if (error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
                    for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                      element_error = std::max(element_error, std::abs(temperr[i]));
                  else if (error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
                    for (unsigned int i=0; i != 6; ++i)
                      element_error = std::max(element_error, std::abs(temperr[i]));
                  else if (error_estimator.error_norm.type(var) == L2)
                    element_error += JxW[sp]*TensorTools::norm_sq(temperr[0]);
                  else if (error_estimator.error_norm.type(var) == H1_SEMINORM)
                    for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                      element_error += JxW[sp]*TensorTools::norm_sq(temperr[i]);
                  else if (error_estimator.error_norm.type(var) == H1_X_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(temperr[0]);
                  else if (error_estimator.error_norm.type(var) == H1_Y_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(temperr[1]);
                  else if (error_estimator.error_norm.type(var) == H1_Z_SEMINORM)
                    element_error += JxW[sp]*TensorTools::norm_sq(temperr[2]);
                  else if (error_estimator.error_norm.type(var) == H2_SEMINORM)
                    {
                      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
                        element_error += JxW[sp]*TensorTools::norm_sq(temperr[i]);
                      // Off diagonal terms enter into the Hessian norm twice
                      for (unsigned int i=3; i != 6; ++i)
                        element_error += JxW[sp]*2*TensorTools::norm_sq(temperr[i]);
                    }

                } // End loop over sample points

              if (error_estimator.error_norm.type(var) == L_INF ||
                  error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                  error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
                new_error_per_cell[e] += error_estimator.error_norm.weight(var) * element_error;
              else if (error_estimator.error_norm.type(var) == L2 ||
                       error_estimator.error_norm.type(var) == H1_SEMINORM ||
                       error_estimator.error_norm.type(var) == H1_X_SEMINORM ||
                       error_estimator.error_norm.type(var) == H1_Y_SEMINORM ||
                       error_estimator.error_norm.type(var) == H1_Z_SEMINORM ||
                       error_estimator.error_norm.type(var) == H2_SEMINORM)
                new_error_per_cell[e] += error_estimator.error_norm.weight_sq(var) * element_error;
              else
                libmesh_error_msg("Unsupported error norm type!");
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

      if(this->error_estimator.patch_reuse)
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
                error_per_cell[e_p_id] =
                  static_cast<ErrorVectorReal>(std::sqrt(new_error_per_cell[i]));
            }
          else
            {
              libmesh_assert (error_estimator.error_norm.type(0) == L_INF ||
                              error_estimator.error_norm.type(0) == W1_INF_SEMINORM ||
                              error_estimator.error_norm.type(0) == W2_INF_SEMINORM);
              Threads::spin_mutex::scoped_lock acquire(Threads::spin_mtx);
              if (!error_per_cell[e_p_id])
                error_per_cell[e_p_id] =
                  static_cast<ErrorVectorReal>(new_error_per_cell[i]);
            }

        } // End loop over every element in patch

    } // end element loop

} // End () operator definition

} // namespace libMesh
