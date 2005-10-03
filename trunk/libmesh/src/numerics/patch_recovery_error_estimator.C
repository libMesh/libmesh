// $Id: patch_recovery_error_estimator.C,v 1.14 2005-10-03 18:55:53 spetersen Exp $

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
#include "quadrature_gauss.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "quadrature_grid.h"
#include "system.h"
#include "mesh.h"



//-----------------------------------------------------------------
// PatchRecoveryErrorEstimator implementations
std::vector<Real> PatchRecoveryErrorEstimator::specpoly(const unsigned int dim,
							const Order order,
							const Real x,
							const Real y,
							const Real z,
							const unsigned int matsize)
{
  std::vector<Real> psi;
  psi.reserve(matsize);
    
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
	      for (int yexp=poly_deg-xexp; yexp >= 0; yexp--) // use an int for yexp since we -- it
		for (int zexp=poly_deg-xexp-yexp; zexp >=0; zexp--) // use an int for zexp since we -- it
		  psi.push_back(std::pow(x,xexp)*std::pow(y,yexp)*std::pow(z,zexp));
	    break;
	  }

	  // 2D spectral polynomial basis functions
	case 2:
	  {
	    for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
	      for (int yexp=poly_deg-xexp; yexp >= 0; yexp--) // use an int for yexp since we -- it
		psi.push_back(std::pow(x,xexp)*std::pow(y,yexp));
	    break;
	  }

	  // 1D spectral polynomial basis functions
	case 1:
	  {
	    for (int xexp=poly_deg; xexp >= 0; xexp--) // use an int for xexp since we -- it
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
						  std::vector<float>& error_per_cell)
{
  START_LOG("estimate_error()", "PatchRecoveryErrorEstimator");

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
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);


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

      const Order element_order (elem->default_order());
      
      // Build a patch containing the current element
      // and its neighbors on the local processor
      std::set<const Elem*> patch;


      this->build_patch_from_local_neighbors (elem, patch);


      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
	{
	  // Possibly skip this variable
	  if (!component_scale.empty())
	    if (component_scale[var] == 0.0) continue;
	  
	  // The type of finite element to use for this variable
	  const FEType& fe_type = dof_map.variable_type (var);

	  // Finite element object for use in this patch
	  AutoPtr<FEBase> fe (FEBase::build (dim, fe_type));
	  
	  // Build an appropriate Gaussian quadrature rule
	  QGauss qrule (dim, fe_type.default_quadrature_order());

	  // Tell the finite element about the quadrature rule.
	  fe->attach_quadrature_rule (&qrule);
      
	  // Get Jacobian values, etc..
	  const std::vector<Real>&                       JxW     = fe->get_JxW();
	  const std::vector<Point>&                      q_point = fe->get_xyz();
	  const std::vector<std::vector<Real> >&         phi     = fe->get_phi();
	  const std::vector<std::vector<RealGradient> >& dphi    = fe->get_dphi();
      
	  // global DOF indices
	  std::vector<unsigned int> dof_indices;

	  unsigned int matsize=1;

	  // Computes the approprite size for the patch projection matrices
	  // and vectors; 
	  for (unsigned int pascal_level=1;
	       pascal_level<static_cast<unsigned int>(element_order)+1;
	       pascal_level++)
	    matsize += this->factorial(pascal_level+dim-1)/this->factorial(pascal_level);
	  	  
	  DenseMatrix<Number> Kp(matsize,matsize);
	  DenseVector<Number>
	    Fx(matsize), Pu_x_h(matsize),
	    Fy(matsize), Pu_y_h(matsize),
	    Fz(matsize), Pu_z_h(matsize);

	  //------------------------------------------------------
	  // Loop over each element in the patch and compute their
	  // contribution to the patch gradient projection.
	  std::set<const Elem*>::const_iterator        patch_it  = patch.begin();
	  const std::set<const Elem*>::const_iterator  patch_end = patch.end();

	  for (; patch_it != patch_end; ++patch_it)
	    {
	      // The pth element in the patch
	      const Elem* e_p = *patch_it;

	      // Reinitialize the finite element data for this element
	      fe->reinit (e_p);

	      // Get the global DOF indices for the current variable
	      // in the current element
	      dof_map.dof_indices (e_p, dof_indices, var);

	      // Huh? Something is horribly WRONG!
	      assert (dof_indices.size() == phi.size());

	      const unsigned int n_dofs = dof_indices.size();
	      const unsigned int n_qp   = qrule.n_points();

	      // Compute the projection components from this cell.
	      // \int_{Omega_e} \psi_i \psi_j = \int_{Omega_e} du_h/dx_k \psi_i
	      for (unsigned int qp=0; qp<n_qp; qp++)
		{
		  // The x,y,z location of the current quadrature point
		  const Real
		    x = q_point[qp](0),
		    y = q_point[qp](1),
		    z = q_point[qp](2);
		    
		  // Compute the gradient on the current patch element
		  // at the quadrature point
		  Gradient grad_u_h;

		  for (unsigned int i=0; i<n_dofs; i++)
		    // grad_u_h += dphi[i][qp]*system.current_solution(dof_indices[i]);
		    grad_u_h.add_scaled (dphi[i][qp],
					 system.current_solution(dof_indices[i]));

		  // Construct the shape function values for the patch projection
		  std::vector<Real> psi(specpoly(dim, element_order, x, y, z, matsize));
		  
		  // Patch matrix contribution
		  for (unsigned int i=0; i<Kp.m(); i++)
		    for (unsigned int j=0; j<Kp.n(); j++)
		      Kp(i,j) += JxW[qp]*psi[i]*psi[j];
		  
		  // Patch RHS contributions
		  for (unsigned int i=0; i<psi.size(); i++)
		    {
		      Fx(i) += JxW[qp]*grad_u_h(0)*psi[i];
		      Fy(i) += JxW[qp]*grad_u_h(1)*psi[i];
		      Fz(i) += JxW[qp]*grad_u_h(2)*psi[i];
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
	  

	  
	  //--------------------------------------------------
	  // Finally, estimate the error in the current variable
	  // for the current element by computing ||Pgrad_u_h - grad_u_h||

	  fe->reinit(elem);
	  //reinitialize element
	  
	  dof_map.dof_indices (elem, dof_indices, var);
	  const unsigned int n_dofs = dof_indices.size();
	  
	  // For linear elments, grad is a constant, so we need to compute
	  // grad u_h once on the element.  Also as G_H u_h - gradu_h is linear
	  // on an element, it assumes its maximum at a vertex of the element
	  Real error = 0;
	  if (element_order == FIRST)
	    {
	      // compute the gradient once at any quadrature point(as its constant)
	      Gradient grad_u_h;
	      for (unsigned int j=0; j<n_dofs; j++)
		grad_u_h.add_scaled (dphi[j][0],
				     system.current_solution(dof_indices[j]));
	      // loop over element vertices
	      for (unsigned int n=0; n<elem->n_nodes(); n++)
		{
		  // Real temperrx=0,temperry=0,temperrz=0;
		  Number temperrx=0,temperry=0,temperrz=0;
		  const Point nodpt=elem->point(n);
		  const Real
		    x = nodpt(0),
		    y = nodpt(1),
		    z = nodpt(2);
		  
		  std::vector<Real> psi(specpoly(dim, element_order,x,y,z,matsize));
		  // get psi-basis values at vertex
		  
		  for (unsigned int i=0; i<matsize; i++)
		    {
		      temperrx += psi[i]*Pu_x_h(i);
		      temperry += psi[i]*Pu_y_h(i);
		      temperrz += psi[i]*Pu_z_h(i);
		    }
		  temperrx -= grad_u_h(0);
		  temperry -= grad_u_h(1);
		  temperrz -= grad_u_h(2);
		  
		  // temperrx = std::abs(temperrx);
		  // temperry = std::abs(temperry);
		  // temperrz = std::abs(temperrz);
		  
		  // error = std::max(temperrz,std::max(temperry,temperrx));

		  error = std::max(std::abs(temperrz),
				   std::max(std::abs(temperry),
					    std::abs(temperrx)));

		} // end vertex loop
	    } // end piecewise linear case
	  else
	    {
	      // we approximate the max norm by sampling over a set of points
	      // in future we may add specialized routines for specific cases
	      // or use some optimization package
	      const Order qorder = TENTH;
	      // build a "fake" quadrature rule for the element
	      QGrid samprule (dim, qorder);
	      fe->attach_quadrature_rule (&samprule);
	      fe->reinit(elem);
	      
	      //const std::vector<Real>&                       JxW2     = fe->get_JxW();
	      const std::vector<Point>&                      samppt   = fe->get_xyz();
	      //const std::vector<std::vector<Real> >&         phi2     = fe->get_phi();
	      const std::vector<std::vector<RealGradient> >& dphi2    = fe->get_dphi();

	      const unsigned int n_sp = samprule.n_points();
	      for (unsigned int sp=0; sp< n_sp; sp++)
		{
		  // Real temperrx=0,temperry=0,temperrz=0;
		  Number temperrx=0,temperry=0,temperrz=0;
		  const Real
		    x = samppt[sp](0),
		    y = samppt[sp](1),
		    z = samppt[sp](2);
		  
		  // Comput the gradient at the current sample point
		  Gradient grad_u_h;
		  
		  for (unsigned int i=0; i<n_dofs; i++)
		    grad_u_h.add_scaled (dphi2[i][sp],
				       system.current_solution(dof_indices[i]));
		  // Compute the phi values at the current sample point
		  std::vector<Real> psi(specpoly(dim,element_order, x,y,z, matsize));
		  for (unsigned int i=0; i<matsize; i++)
		    {
		      temperrx += psi[i]*Pu_x_h(i);
		      temperry += psi[i]*Pu_y_h(i);
		      temperrz += psi[i]*Pu_z_h(i);
		    }
		  temperrx -= grad_u_h(0);
		  temperry -= grad_u_h(1);
		  temperrz -= grad_u_h(2);

		  // temperrx = std::abs(temperrx);
		  // temperry = std::abs(temperry);
		  // temperrz = std::abs(temperrz);
		  
		  // error = std::max(temperrz,std::max(temperry,temperrx));		  

		  error = std::max(std::abs(temperrz),
				   std::max(std::abs(temperry),
					    std::abs(temperrx)));

		} // end sample_point_loop
	    } // end P>1 Error loop
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



void PatchRecoveryErrorEstimator::build_patch_from_local_neighbors (const Elem* e0,
								    std::set<const Elem*> patch,
								    const unsigned int target_patch_size)
{
  START_LOG("build_patch_from_local_neighbors()", "PatchRecoveryErrorEstimator");
  
  // Make sure we are building a patch for an active, local element.
  // (Are these restrictions necessary?)
  assert (e0 != NULL);
  assert (e0->processor_id() == libMesh::processor_id());
  assert (e0->active());
  
  // First add the element of interest.
  patch.insert (e0);

  // Repeatedly add the neighbors of the elements in the patch until
  // the target patch size is met
  while (patch.size() < target_patch_size)
    {
      // It is possible that the target patch size is larger than the number
      // of active elements in the mesh.  Since we don't have access to the
      // Mesh object here the only way we can detect this case is by detecting
      // a "stagnant patch," i.e. a patch whose size does not increase after adding
      // face neighbors
      const unsigned int old_patch_size = patch.size();
      
      // Loop over all the elements in the patch
      std::set<const Elem*>::const_iterator       it  = patch.begin();
      const std::set<const Elem*>::const_iterator end = patch.end();

      for (; (it != end) && (patch.size() < target_patch_size); ++it)
	{
	  // Convenience.  Keep the syntax simple.
	  const Elem* elem = *it;

	  for (unsigned int s=0; s<elem->n_sides(); s++)
	    if (elem->neighbor(s) != NULL)        // we have a neighbor on this side
	      {
		const Elem* neighbor = elem->neighbor(s);
		
		if (neighbor->active())           // ... and that neighbor is active
		  if (neighbor->processor_id() ==
		      libMesh::processor_id())    // ... and belongs to this processor
		    patch.insert (neighbor);      // ... then add it to the patch
		  
		else                              // ... the neighbor is *not* active,
		  {                               // ... so add *all* its active, local children to the patch
		    std::vector<const Elem*> active_children;

		    neighbor->active_family_tree (active_children);

		    for (unsigned int c=0; c<active_children.size(); c++)
		      if (active_children[c]->processor_id() == libMesh::processor_id())
			patch.insert (active_children[c]);
		  }
	      }
	}
      
      // Check for a "stagnant" patch
      if (patch.size() == old_patch_size)
	{
	  std::cerr << "ERROR: stagnant patch of "
		    << patch.size() << " elements."
		    << std::endl
		    << "Does your target patch size exceed the number of elements in the mesh?"
		    << std::endl;
	  here();
	  break;
	}
    } // end while loop

  
  // make sure all the elements in the patch are active and local
  // if we are in debug mode
#ifdef DEBUG
  {
    std::set<const Elem*>::const_iterator       it  = patch.begin();
    const std::set<const Elem*>::const_iterator end = patch.end();
    
    for (; it != end; ++it)
      {
	// Convenience.  Keep the syntax simple.
	const Elem* elem = *it;

	assert (elem->active());
	assert (elem->processor_id() == libMesh::processor_id());
      }
  }
#endif

  STOP_LOG("build_patch_from_local_neighbors()", "PatchRecoveryErrorEstimator");
}
