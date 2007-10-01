// $Id: patch_recovery_error_estimator.C,v 1.27 2007-10-01 23:17:06 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
		  psi.push_back(std::pow(x,xexp)*std::pow(y,yexp)*std::pow(z,zexp));
                }
	    break;
	  }

	  // 2D spectral polynomial basis functions
	case 2:
	  {
	    for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
              {
                int yexp = poly_deg - xexp;
		psi.push_back(std::pow(x,xexp)*std::pow(y,yexp));
              }
	    break;
	  }

	  // 1D spectral polynomial basis functions
	case 1:
	  {
            int xexp = poly_deg;
	    psi.push_back(std::pow(x,xexp));
	    break;
	  }
	  
	default:
	  error();
	}
    }

  return psi;
}
    
  

void PatchRecoveryErrorEstimator::estimate_error (const System& system,
						  ErrorVector& error_per_cell,
					          bool)
{
  START_LOG("estimate_error()", "PatchRecoveryErrorEstimator");

#ifdef ENABLE_SECOND_DERIVATIVES
  assert (_sobolev_order == 1 || _sobolev_order == 2);
#else
  assert (_sobolev_order == 1);
#endif

  // The current mesh
  const Mesh& mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();
  
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.max_elem_id());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Check for the use of component_mask
  this->convert_component_mask_to_scale();

  // Check for a valid component_scale
  if (!component_scale.empty())
    if (component_scale.size() != n_vars)
      {
	std::cerr << "ERROR: component_scale is the wrong size:"
		  << std::endl
		  << " component_scale.size()=" << component_scale.size()
		  << std::endl
		  << ", n_vars=" << n_vars
		  << std::endl;
	error();
      }


  //------------------------------------------------------------
  // Iterate over all the active elements in the mesh
  // that live on this processor.
  MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 
  
  for (; elem_it != elem_end; ++elem_it)
    {
      // elem is necessarily an active element on the local processor
      const Elem* elem = *elem_it;

      // Build a patch containing the current element
      // and its neighbors on the local processor
      Patch patch;

      // We log the time spent building patches separately
      PAUSE_LOG("estimate_error()", "PatchRecoveryErrorEstimator");

      // Use user specified patch size and growth strategy
      patch.build_around_element (elem, target_patch_size,
				  patch_growth_strategy);

      RESTART_LOG("estimate_error()", "PatchRecoveryErrorEstimator");

      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
	{
	  // Possibly skip this variable
	  if (!component_scale.empty())
	    if (component_scale[var] == 0.0) continue;
	  
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
	  const std::vector<std::vector<Real> >&         phi     = fe->get_phi();
	  const std::vector<std::vector<RealGradient> > *dphi = NULL;
          if (_sobolev_order == 1)
            dphi = &(fe->get_dphi());
#ifdef ENABLE_SECOND_DERIVATIVES
	  const std::vector<std::vector<RealTensor> >  *d2phi = NULL;
          if (_sobolev_order == 2)
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
	  DenseVector<Number>
	    Fx(matsize), Pu_x_h(matsize), // Also xx
	    Fy(matsize), Pu_y_h(matsize), // Also yy
	    Fz(matsize), Pu_z_h(matsize); // Also zz
          DenseVector<Number> Fxy, Pu_xy_h, Fxz, Pu_xz_h, Fyz, Pu_yz_h;
          if (_sobolev_order == 2)
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
	      assert (dof_indices.size() == phi.size());

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

                  if (_sobolev_order == 1)
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
                  else if (_sobolev_order == 2)
                    {
#ifdef ENABLE_SECOND_DERIVATIVES
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
		      std::cerr << "ERROR:  --enable-second-derivatives is required\n"
				<< "        for _sobolev_order == 2!\n";
		      error();
#endif
                    }
		} // end quadrature loop
	    } // end patch loop
	  

	  
	  //--------------------------------------------------
	  // Now we have fully assembled the projection system
	  // for this patch.  Project the gradient components.
	  // MAY NEED TO USE PARTIAL PIVOTING!
	  Kp.lu_solve (Fx, Pu_x_h);
	  Kp.lu_solve (Fy, Pu_y_h);
	  Kp.lu_solve (Fz, Pu_z_h);
          if (_sobolev_order == 2)
            {
              Kp.lu_solve(Fxy, Pu_xy_h);
              Kp.lu_solve(Fxz, Pu_xz_h);
              Kp.lu_solve(Fyz, Pu_yz_h);
            }
	  
	  //--------------------------------------------------
	  // Finally, estimate the error in the current variable
	  // for the current element by computing ||P grad_u_h - grad_u_h||
          // or ||P hess_u_h - hess_u_h|| in the infinity (max) norm

	  fe->reinit(elem);
	  //reinitialize element
	  
	  dof_map.dof_indices (elem, dof_indices, var);
	  const unsigned int n_dofs = dof_indices.size();
	  
	  // For linear elments, grad is a constant, so we need to compute
	  // grad u_h once on the element.  Also as G_H u_h - gradu_h is linear
	  // on an element, it assumes its maximum at a vertex of the element
	  Real error = 0;
	  // we approximate the max norm by sampling over a set of points
	  // in future we may add specialized routines for specific cases
	  // or use some optimization package
	  const Order qorder = element_order;

	  // build a "fake" quadrature rule for the element
	  QGrid samprule (dim, qorder);
	  fe->attach_quadrature_rule (&samprule);
	  fe->reinit(elem);
	      
	  const unsigned int n_sp = samprule.n_points();
	  for (unsigned int sp=0; sp< n_sp; sp++)
	    {
	      std::vector<Number> temperr(6,0.0); // x,y,z or xx,yy,zz,xy,xz,yz
	  
              if (_sobolev_order == 1)
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
              else if (_sobolev_order == 2)
                {
#ifdef ENABLE_SECOND_DERIVATIVES
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
		      std::cerr << "ERROR:  --enable-second-derivatives is required\n"
				<< "        for _sobolev_order == 2!\n";
		      error();
#endif
                }
              for (unsigned int i=0; i != 6; ++i)
                error = std::max(error, std::abs(temperr[i]));

	    } // end sample_point_loop
	  const int e_id=elem->id();
	  error_per_cell[e_id] += error;	  
	} // end variable loop  
    } // end element loop

  // Each processor has now computed the error contribuions
  // for its local elements, and error_per_cell contains 0 for all the
  // non-local elements.  Summing the vector will provide the L_oo value
  // for each element, local or remote
  this->reduce_error(error_per_cell);
  
  STOP_LOG("estimate_error()", "PatchRecoveryErrorEstimator");
}
