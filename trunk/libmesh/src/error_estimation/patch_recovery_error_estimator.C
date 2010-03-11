// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <cmath>     // for std::sqrt std::pow std::abs


// Local Includes
#include "libmesh_common.h"
#include "patch_recovery_error_estimator.h"
#include "dof_map.h"
#include "fe.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "error_vector.h"
#include "quadrature_gauss.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "patch.h"
#include "quadrature_grid.h"
#include "system.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "tensor_value.h"


//-----------------------------------------------------------------
// PatchRecoveryErrorEstimator implementations
std::vector<Real> PatchRecoveryErrorEstimator::specpoly(const unsigned int dim,
							const Order order,
							const Point p,
							const unsigned int matsize)
{
  std::vector<Real> psi;
  psi.reserve(matsize);
  Real x = p(0), y=0., z=0.;
  if (dim > 1)
    y = p(1);
  if (dim > 2)
    z = p(2);
    
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
		  psi.push_back(std::pow(x,static_cast<Real>(xexp))*
				std::pow(y,static_cast<Real>(yexp))*
				std::pow(z,static_cast<Real>(zexp)));
                }
	    break;
	  }

	  // 2D spectral polynomial basis functions
	case 2:
	  {
	    for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
              {
                int yexp = poly_deg - xexp;
		psi.push_back(std::pow(x,static_cast<Real>(xexp))*
			      std::pow(y,static_cast<Real>(yexp)));
              }
	    break;
	  }

	  // 1D spectral polynomial basis functions
	case 1:
	  {
            int xexp = poly_deg;
	    psi.push_back(std::pow(x,static_cast<Real>(xexp)));
	    break;
	  }
	  
	default:
	  libmesh_error();
	}
    }

  return psi;
}
    
  

void PatchRecoveryErrorEstimator::estimate_error (const System& system,
						  ErrorVector& error_per_cell,
					          const NumericVector<Number>* solution_vector,
					          bool)
{
  START_LOG("estimate_error()", "PatchRecoveryErrorEstimator");

  // The current mesh
  const MeshBase& mesh = system.get_mesh();
  
  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.max_elem_id());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Prepare current_local_solution to localize a non-standard
  // solution vector if necessary
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number>* newsol =
        const_cast<NumericVector<Number>*>(solution_vector);
      System &sys = const_cast<System&>(system);
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

  // Each processor has now computed the error contribuions
  // for its local elements, and error_per_cell contains 0 for all the
  // non-local elements.  Summing the vector will provide the L_oo value
  // for each element, local or remote
  this->reduce_error(error_per_cell);
  
  // If we used a non-standard solution before, now is the time to fix
  // the current_local_solution
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number>* newsol =
        const_cast<NumericVector<Number>*>(solution_vector);
      System &sys = const_cast<System&>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }

  STOP_LOG("estimate_error()", "PatchRecoveryErrorEstimator");
}



void PatchRecoveryErrorEstimator::EstimateError::operator()(const ConstElemRange &range) const
{
  // The current mesh
  const MeshBase& mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();
  
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();


  //------------------------------------------------------------
  // Iterate over all the elements in the range.
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it!=range.end(); ++elem_it)
    {
      // elem is necessarily an active element on the local processor
      const Elem* elem = *elem_it;

      // We'll need an index into the error vector
      const int e_id=elem->id();

      // Build a patch containing the current element
      // and its neighbors on the local processor
      Patch patch;

      // Use user specified patch size and growth strategy
      patch.build_around_element (elem, error_estimator.target_patch_size,
				  error_estimator.patch_growth_strategy);

      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
	{
#ifndef DEBUG
  #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          libmesh_assert (error_estimator.error_norm.type(var) == L2 ||
			  error_estimator.error_norm.type(var) == H1_SEMINORM ||
                          error_estimator.error_norm.type(var) == H2_SEMINORM ||
                          error_estimator.error_norm.type(var) == L_INF ||
                          error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                          error_estimator.error_norm.type(var) == W2_INF_SEMINORM);
  #else
          libmesh_assert (error_estimator.error_norm.type(var) == L2 ||
			  error_estimator.error_norm.type(var) == L_INF ||
			  error_estimator.error_norm.type(var) == H1_SEMINORM ||
                          error_estimator.error_norm.type(var) == W1_INF_SEMINORM);
  #endif
          if (var > 0)
            // We can't mix L_inf and L_2 norms
            libmesh_assert (((error_estimator.error_norm.type(var) == L2 ||
			      error_estimator.error_norm.type(var) == H1_SEMINORM ||
                              error_estimator.error_norm.type(var) == H2_SEMINORM) &&
                             (error_estimator.error_norm.type(var-1) == L2 ||
			      error_estimator.error_norm.type(var-1) == H1_SEMINORM ||
                              error_estimator.error_norm.type(var-1) == H2_SEMINORM)) ||
                            ((error_estimator.error_norm.type(var) == L_INF ||
                              error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
                              error_estimator.error_norm.type(var) == W2_INF_SEMINORM) &&
                             (error_estimator.error_norm.type(var-1) == L_INF ||
                              error_estimator.error_norm.type(var-1) == W1_INF_SEMINORM ||
                              error_estimator.error_norm.type(var-1) == W2_INF_SEMINORM)));
#endif

	  // Possibly skip this variable
	  if (error_estimator.error_norm.weight(var) == 0.0) continue;
	  
	  // The type of finite element to use for this variable
	  const FEType& fe_type = dof_map.variable_type (var);

          const Order element_order  = static_cast<Order>
            (fe_type.order + elem->p_level());
      
	  // Finite element object for use in this patch
	  AutoPtr<FEBase> fe (FEBase::build (dim, fe_type));
	  
	  // Build an appropriate Gaussian quadrature rule
	  AutoPtr<QBase> qrule (fe_type.default_quadrature_rule(dim));

	  // Tell the finite element about the quadrature rule.
	  fe->attach_quadrature_rule (qrule.get());
      
	  // Get Jacobian values, etc..
	  const std::vector<Real>&                       JxW     = fe->get_JxW();
	  const std::vector<Point>&                      q_point = fe->get_xyz();

          // Get whatever phi/dphi/d2phi values we need.  Avoid
          // getting them unless the requested norm is actually going
          // to use them.

	  const std::vector<std::vector<Real> >         *phi = NULL;
          // If we're using phi to assert the correct dof_indices
          // vector size later, then we'll need to get_phi whether we
          // plan to use it or not.
#ifdef NDEBUG
          if (error_estimator.error_norm.type(var) == L2 ||
              error_estimator.error_norm.type(var) == L_INF)
#endif
            phi = &(fe->get_phi());

	  const std::vector<std::vector<RealGradient> > *dphi = NULL;
          if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
              error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
            dphi = &(fe->get_dphi());

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	  const std::vector<std::vector<RealTensor> >  *d2phi = NULL;
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            d2phi = &(fe->get_d2phi());
#endif
      
	  // global DOF indices
	  std::vector<unsigned int> dof_indices;

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
          if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
              error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
              error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
              Fx.resize(matsize); Pu_x_h.resize(matsize); // stores xx in W2 cases
              Fy.resize(matsize); Pu_y_h.resize(matsize); // stores yy in W2 cases
              Fz.resize(matsize); Pu_z_h.resize(matsize); // stores zz in W2 cases
            }
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
              Fxy.resize(matsize); Pu_xy_h.resize(matsize);
              Fxz.resize(matsize); Pu_xz_h.resize(matsize);
              Fyz.resize(matsize); Pu_yz_h.resize(matsize);
            }

	  //------------------------------------------------------
	  // Loop over each element in the patch and compute their
	  // contribution to the patch gradient projection.
	  Patch::const_iterator        patch_it  = patch.begin();
	  const Patch::const_iterator  patch_end = patch.end();

	  for (; patch_it != patch_end; ++patch_it)
	    {
	      // The pth element in the patch
	      const Elem* e_p = *patch_it;

	      // Reinitialize the finite element data for this element
	      fe->reinit (e_p);

	      // Get the global DOF indices for the current variable
	      // in the current element
	      dof_map.dof_indices (e_p, dof_indices, var);
	      libmesh_assert (dof_indices.size() == phi->size());

	      const unsigned int n_dofs = dof_indices.size();
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
		      for (unsigned int i=0; i<psi.size(); i++)			
			  F(i) = JxW[qp]*u_h*psi[i];			

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
		      for (unsigned int i=0; i<psi.size(); i++)
		        {
		          Fx(i) += JxW[qp]*grad_u_h(0)*psi[i];
		          Fy(i) += JxW[qp]*grad_u_h(1)*psi[i];
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
		      for (unsigned int i=0; i<psi.size(); i++)
		        {
		          Fx(i)  += JxW[qp]*hess_u_h(0,0)*psi[i];
		          Fy(i)  += JxW[qp]*hess_u_h(1,1)*psi[i];
		          Fz(i)  += JxW[qp]*hess_u_h(2,2)*psi[i];
		          Fxy(i) += JxW[qp]*hess_u_h(0,1)*psi[i];
		          Fxz(i) += JxW[qp]*hess_u_h(0,2)*psi[i];
		          Fyz(i) += JxW[qp]*hess_u_h(1,2)*psi[i];
		        }
#else
		      libMesh::err << "ERROR:  --enable-second-derivatives is required\n"
				    << "        for _sobolev_order == 2!\n";
		      libmesh_error();
#endif
                    }
                  else
	            libmesh_error();
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
          if (error_estimator.error_norm.type(var) == H1_SEMINORM ||
              error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
              error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
	      Kp.lu_solve (Fx, Pu_x_h);
	      Kp.lu_solve (Fy, Pu_y_h);
	      Kp.lu_solve (Fz, Pu_z_h);
            }
          if (error_estimator.error_norm.type(var) == H2_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
            {
              Kp.lu_solve(Fxy, Pu_xy_h);
              Kp.lu_solve(Fxz, Pu_xz_h);
              Kp.lu_solve(Fyz, Pu_yz_h);
            }
	  
	  //--------------------------------------------------
	  // Finally, estimate the error in the current variable
	  // for the current element by computing ||P u_h - u_h|| or ||P grad_u_h - grad_u_h||
	  // or ||P hess_u_h - hess_u_h|| according to the requested
	  // seminorm

	  dof_map.dof_indices (elem, dof_indices, var);
	  const unsigned int n_dofs = dof_indices.size();
	  
	  Real element_error = 0;

	  // we approximate the max norm by sampling over a set of points
	  // in future we may add specialized routines for specific cases
	  // or use some optimization package
	  const Order qorder = element_order;

	  // build a "fake" quadrature rule for the element
	  //
	  // For linear elements, grad is a constant, so we need
	  // to compute grad u_h once on the element.  Also as G_H
	  // u_h - gradu_h is linear on an element, it assumes its
	  // maximum at a vertex of the element

	  QGrid samprule (dim, qorder);

          if (error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
	    fe->attach_quadrature_rule (&samprule);

	  // reinitialize element for integration or sampling
	  fe->reinit(elem);
	      
	  const unsigned int n_sp = JxW.size();
	  for (unsigned int sp=0; sp< n_sp; sp++)
	    {
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
	              temperr[1] += psi[i]*Pu_y_h(i);
	              temperr[2] += psi[i]*Pu_z_h(i);
	            }
	          temperr[0] -= grad_u_h(0);
	          temperr[1] -= grad_u_h(1);
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
	              temperr[1] += psi[i]*Pu_y_h(i);
	              temperr[2] += psi[i]*Pu_z_h(i);
	              temperr[3] += psi[i]*Pu_xy_h(i);
	              temperr[4] += psi[i]*Pu_xz_h(i);
	              temperr[5] += psi[i]*Pu_yz_h(i);
	            }
	          temperr[0] -= hess_u_h(0,0);
	          temperr[1] -= hess_u_h(1,1);
	          temperr[2] -= hess_u_h(2,2);
	          temperr[3] -= hess_u_h(0,1);
	          temperr[4] -= hess_u_h(0,2);
	          temperr[5] -= hess_u_h(1,2);
#else
		      libMesh::err << "ERROR:  --enable-second-derivatives is required\n"
				    << "        for _sobolev_order == 2!\n";
		      libmesh_error();
#endif
                }
              if (error_estimator.error_norm.type(var) == L_INF)
                element_error = std::max(element_error, std::abs(temperr[0]));
              else if (error_estimator.error_norm.type(var) == W1_INF_SEMINORM)
                for (unsigned int i=0; i != 3; ++i)
                  element_error = std::max(element_error, std::abs(temperr[i]));
              else if (error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
                for (unsigned int i=0; i != 6; ++i)
                  element_error = std::max(element_error, std::abs(temperr[i]));
              else if (error_estimator.error_norm.type(var) == L2)
                element_error += JxW[sp]*libmesh_norm(temperr[0]);
              else if (error_estimator.error_norm.type(var) == H1_SEMINORM)
                for (unsigned int i=0; i != 3; ++i)
                  element_error += JxW[sp]*libmesh_norm(temperr[i]);
              else if (error_estimator.error_norm.type(var) == H2_SEMINORM)
                {
                  for (unsigned int i=0; i != 3; ++i)
                    element_error += JxW[sp]*libmesh_norm(temperr[i]);
                // Off diagonal terms enter into the Hessian norm twice
                  for (unsigned int i=3; i != 6; ++i)
                    element_error += JxW[sp]*2*libmesh_norm(temperr[i]);
                }

	    } // end sample_point_loop

	  // The patch error estimator works element-by-element --
	  // there is no need to get a mutex on error_per_cell!
	  if (error_estimator.error_norm.type(var) == L_INF ||
	      error_estimator.error_norm.type(var) == W1_INF_SEMINORM ||
              error_estimator.error_norm.type(var) == W2_INF_SEMINORM)
	    error_per_cell[e_id] += error_estimator.error_norm.weight(var) * element_error;	  
          else if (error_estimator.error_norm.type(var) == L2 ||
		   error_estimator.error_norm.type(var) == H1_SEMINORM ||
                   error_estimator.error_norm.type(var) == H2_SEMINORM)
	    error_per_cell[e_id] += error_estimator.error_norm.weight_sq(var) * element_error;	  
          else
	    libmesh_error();
	} // end variable loop  

      if (error_estimator.error_norm.type(0) == L2 ||
	  error_estimator.error_norm.type(0) == H1_SEMINORM ||
          error_estimator.error_norm.type(0) == H2_SEMINORM)
        error_per_cell[e_id] = std::sqrt(error_per_cell[e_id]);

    } // end element loop
}
