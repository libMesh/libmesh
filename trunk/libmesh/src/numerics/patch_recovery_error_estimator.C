// $Id: patch_recovery_error_estimator.C,v 1.3 2004-06-08 14:45:50 spetersen Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>    // for sqrt


// Local Includes
#include "libmesh_common.h"
#include "patch_recovery_error_estimator.h"
#include "dof_map.h"
#include "fe.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "quadrature_gauss.h"
#include "libmesh_logging.h"




//-----------------------------------------------------------------
// ErrorEstimator implementations
void PatchRecoveryErrorEstimator::estimate_error (const SteadySystem& system,
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


  // Check for a valid component_mask
  if (!component_mask.empty())
    if (component_mask.size() != n_vars)
      {
	std::cerr << "ERROR: component_mask is the wrong size:"
		  << std::endl
		  << " component_mask.size()=" << component_mask.size()
		  << std::endl
		  << ", n_vars=" << n_vars
		  << std::endl;
	error();
      }


  //------------------------------------------------------------
  // Iterate over all the active elements in the mesh
  // that live on this processor.
  const_active_local_elem_iterator       elem_it (mesh.elements_begin());
  const const_active_local_elem_iterator elem_end(mesh.elements_end());
  
  for (; elem_it != elem_end; ++elem_it)
    {
      // elem is necessarily an active element on the local processor
      const Elem* elem = *elem_it;

      // Build a patch containing the current element
      // and its neighbors on the local processor
      std::set<const Elem*> patch;

      this->build_patch_from_local_neighbors (elem, patch);


      //------------------------------------------------------------
      // Process each variable in the system using the current patch
      for (unsigned int var=0; var<n_vars; var++)
	{
	  // Possibly skip this variable
	  if (!component_mask.empty())
	    if (component_mask[var] == false) continue;
	  
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

	  // Dense matrix and vectors for patch projection.
	  // THE SIZES OF THESE NEED TO GET SMARTER!
	  DenseMatrix<Number> Kp(dim+1,dim+1);
	  DenseVector<Number>
	    Fx(dim+1), Pu_x_h(dim+1),
	    Fy(dim+1), Pu_y_h(dim+1),
	    Fz(dim+1), Pu_z_h(dim+1);

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
		  std::vector<Real> psi;
		  {
		    psi.reserve (dim+1);
		    
		    psi.push_back(1);
		    psi.push_back(x);
		    psi.push_back(y);
		    
		    if (dim == 3)
		      psi.push_back(z);
		  }
		  
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
	  error();
	  
	} // end variable loop  
    } // end element loop

  
  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector
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
		    << "Does your target patch size exceed the number of elements?"
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
