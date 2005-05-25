// $Id: kelly_error_estimator.C,v 1.11 2005-05-25 17:54:24 jwpeterson Exp $

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
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh_common.h"
#include "kelly_error_estimator.h"
#include "dof_map.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "libmesh_logging.h"
#include "elem.h"



//-----------------------------------------------------------------
// ErrorEstimator implementations
void KellyErrorEstimator::estimate_error (const System& system,
					  std::vector<float>& error_per_cell)
{
  //  START_LOG("flux_jumps()", "KellyErrorEstimator");
  
  /*

  Conventions for assigning the direction of the normal:
  
  - e & f are global element ids
  
  Case (1.) Elements are at the same level, e<f
            Compute the flux jump on the face and
	    add it as a contribution to error_per_cell[e]
	    and error_per_cell[f]
  
                   ----------------------
		  |           |          |
		  |           |    f     |
		  |           |          |
		  |    e      |---> n    | 
		  |           |          |
		  |           |          |
                   ----------------------


   Case (2.) The neighbor is at a higher level.
             Compute the flux jump on e's face and
	     add it as a contribution to error_per_cell[e]
	     and error_per_cell[f]

                   ----------------------
		  |     |     |          |
		  |     |  e  |---> n    |
		  |     |     |          |
		  |-----------|    f     | 
		  |     |     |          |
		  |     |     |          |
		  |     |     |          |
                   ----------------------
  */
   
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
		  << " n_vars=" << n_vars
		  << std::endl;
	error();
      }
  

  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_mask.empty())
	if (component_mask[var] == false) continue;

      // The (string) name of this variable
      const std::string& var_name = system.variable_name(var);
      
      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);
      
      // Finite element objects for the same face from
      // different sides
      AutoPtr<FEBase> fe_e (FEBase::build (dim, fe_type));
      AutoPtr<FEBase> fe_f (FEBase::build (dim, fe_type));

      // Build an appropriate Gaussian quadrature rule
      QGauss qrule (dim-1, fe_type.default_quadrature_order());

      // Tell the finite element for element e about the quadrature
      // rule.  The finite element for element f need not know about it
      fe_e->attach_quadrature_rule (&qrule);
      
      // By convention we will always do the integration
      // on the face of element e.  Get its Jacobian values, etc..
      const std::vector<Real>&  JxW_face     = fe_e->get_JxW();
      const std::vector<Point>& qface_point  = fe_e->get_xyz();
      const std::vector<Point>& face_normals = fe_e->get_normals();

      // The quadrature points on element f.  These will be computed
      // from the quadrature points on element e.
      std::vector<Point> qp_f;
      
      // The shape function gradients on elements e & f
      const std::vector<std::vector<RealGradient> > & dphi_e = fe_e->get_dphi();
      const std::vector<std::vector<RealGradient> > & dphi_f = fe_f->get_dphi();
      
      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices_e;
      std::vector<unsigned int> dof_indices_f;


      
      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* e = *elem_it;
	  const unsigned int e_id = e->id();
	  
	  // Loop over the neighbors of element e
	  for (unsigned int n_e=0; n_e<e->n_neighbors(); n_e++)
	    {
	      if (e->neighbor(n_e) != NULL) // e is not on the boundary
		{
		  const Elem* f           = e->neighbor(n_e);
		  const unsigned int f_id = f->id();
		
		  if (   //-------------------------------------
		      ((f->active()) &&
		       (f->level() == e->level()) &&
		       (e_id < f_id))                 // Case 1.
		    
		      || //-------------------------------------
		    
		      (f->level() < e->level())       // Case 2.
		    
		      )  //-------------------------------------
		    {		    
		      // Update the shape functions on side s_e of
		      // element e
		      START_LOG("fe_e->reinit()", "KellyErrorEstimator");
		      fe_e->reinit (e, n_e);
		      STOP_LOG("fe_e->reinit()", "KellyErrorEstimator");
		      
		      // Build the side
		      START_LOG("construct side", "KellyErrorEstimator");
		      AutoPtr<Elem> side (e->side(n_e));
		      STOP_LOG("construct side", "KellyErrorEstimator");
		      
		      // Get the maximum h for this side
		      START_LOG("side->hmax()", "KellyErrorEstimator");
		      const Real h = side->hmax();
		      STOP_LOG("side->hmax()", "KellyErrorEstimator");
		      
		      // Get the DOF indices for the two elements
		      START_LOG("dof_indices()", "KellyErrorEstimator");
		      dof_map.dof_indices (e, dof_indices_e, var);
		      dof_map.dof_indices (f, dof_indices_f, var);
		      STOP_LOG("dof_indices()", "KellyErrorEstimator");

		      // The number of DOFS on each element
		      const unsigned int n_dofs_e = dof_indices_e.size();
		      const unsigned int n_dofs_f = dof_indices_f.size();

		      // The number of quadrature points
		      const unsigned int n_qp = qrule.n_points();

		      // Find the location of the quadrature points
		      // on element f
		      START_LOG("inverse_map()", "KellyErrorEstimator");
		      FEInterface::inverse_map (dim, fe_type, f, qface_point, qp_f);
		      STOP_LOG("inverse_map()", "KellyErrorEstimator");
		      
		      // Compute the shape functions on element f
		      // at the quadrature points of element e
		      START_LOG("fe_f->reinit()", "KellyErrorEstimator");
		      fe_f->reinit (f, &qp_f);
		      STOP_LOG("fe_f->reinit()", "KellyErrorEstimator");
		      
		      // The error contribution from this face
		      Real error = 1.e-30;

		    
		      START_LOG("jump integral", "KellyErrorEstimator");
		      // loop over the integration points on the face
		      for (unsigned int qp=0; qp<n_qp; qp++)
			{
			  // The solution gradient from each element
			  Gradient grad_e, grad_f;
			
			  // Compute the solution gradient on element e
			  for (unsigned int i=0; i<n_dofs_e; i++)
			    grad_e.add_scaled (dphi_e[i][qp],
					       system.current_solution(dof_indices_e[i]));
			
			  // Compute the solution gradient on element f
			  for (unsigned int i=0; i<n_dofs_f; i++)
			    grad_f.add_scaled (dphi_f[i][qp],
					       system.current_solution(dof_indices_f[i]));
			

			  // The flux jump at the face 
			  const Number jump = (grad_e - grad_f)*face_normals[qp];

			  // The flux jump squared.  If using complex numbers,
			  // std::norm(z) returns |z|^2, where |z| is the modulus of z.
#ifndef USE_COMPLEX_NUMBERS
			  const Real jump2 = jump*jump;
#else
			  const Real jump2 = std::norm(jump);
#endif

			  // Integrate the error on the face.  The error is
			  // scaled by an additional power of h, where h is
			  // the maximum side length for the element.  This
			  // arises in the definition of the indicator.
			  error += JxW_face[qp]*h*jump2;			
			
			} // End quadrature point loop
		      STOP_LOG("jump integral", "KellyErrorEstimator");
		      
		      // Add the error contribution to elements e & f
                      assert(e_id < error_per_cell.size());
                      assert(f_id < error_per_cell.size());
		      error_per_cell[e_id] += error;
		      error_per_cell[f_id] += error;
		    } // end if case 1 or case 2
		} // if (e->neigbor(n_e) != NULL)

	      // Otherwise, e is on the boundary.  If it happens to
	      // be on a Dirichlet boundary, we need not do anything.
	      // On the other hand, if e is on a Neumann (flux) boundary
	      // with grad(u).n = g, we need to compute the additional residual
	      // (h * \int |g - grad(u_h).n|^2 dS)^(1/2).
	      // We can only do this with some knowledge of the boundary
	      // conditions, i.e. the user must have attached an appropriate
	      // BC function.
	      else
		{
		  if (this->_bc_function != NULL)
		    {
		      START_LOG("boundary integrals", "KellyErrorEstimator");
		      // here();
		  
		      // Update the shape functions on side s_e of element e
		      fe_e->reinit (e, n_e);

		      // The reinitialization also recomputes the locations of
		      // the quadrature points on the side.  By checking if the
		      // first quadrature point on the side is on a flux boundary
		      // for a particular variable, we will determine if the whole
		      // element is on a flux boundary (assuming quadrature points
		      // are strictly contained in the side).
		      if (this->_bc_function(system, qface_point[0], var_name).first)
			{
			  // Build the side
			  AutoPtr<Elem> side (e->side(n_e));

			  // Get the maximum h for this side
			  const Real h = side->hmax();
		    
			  // Get the DOF indices 
			  dof_map.dof_indices (e, dof_indices_e, var);

			  // The number of DOFS on each element
			  const unsigned int n_dofs_e = dof_indices_e.size();

			  // The number of quadrature points
			  const unsigned int n_qp = qrule.n_points();

			  // The error contribution from this face
			  Real error = 1.e-10;
		    
			  // loop over the integration points on the face.
			  for (unsigned int qp=0; qp<n_qp; qp++)
			    {
			      // Value of the imposed flux BC at this quadrature point.
			      const std::pair<bool,Real> flux_bc =
				this->_bc_function(system, qface_point[qp], var_name);

			      // Be sure the BC function still thinks we're on the 
			      // flux boundary.
			      assert (flux_bc.first == true);
			      
			      // The solution gradient from each element
			      Gradient grad_e;
			
			      // Compute the solution gradient on element e
			      for (unsigned int i=0; i<n_dofs_e; i++)
				grad_e.add_scaled (dphi_e[i][qp],
						   system.current_solution(dof_indices_e[i]));

			      // The difference between the desired BC and the approximate solution. 
			      const Number jump = flux_bc.second - grad_e*face_normals[qp];

			      // The flux jump squared.  If using complex numbers,
			      // std::norm(z) returns |z|^2, where |z| is the modulus of z.
#ifndef USE_COMPLEX_NUMBERS
			      const Real jump2 = jump*jump;
#else
			      const Real jump2 = std::norm(jump);
#endif

			      
// 			      std::cout << "Error contribution from "
// 					<< var_name
// 					<< " flux BC: "
// 					<< JxW_face[qp]*h*jump2
// 					<< std::endl;

			      
			      // Integrate the error on the face.  The error is
			      // scaled by an additional power of h, where h is
			      // the maximum side length for the element.  This
			      // arises in the definition of the indicator.
			      error += JxW_face[qp]*h*jump2;			
			
			    } // End quadrature point loop

			  // Add the error contribution to elements e & f
                          assert(e_id < error_per_cell.size());
			  error_per_cell[e_id] += error;
			  
			} // end if side on flux boundary
		      
		      STOP_LOG("boundary integrals", "KellyErrorEstimator");
		    } // end if _bc_function != NULL
		} // end if (e->neighbor(n_e) == NULL)
	    } // end loop over neighbors
	} // End loop over active local elements
    } // End loop over variables


  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector
  // and then take the square-root of each component.  Note
  // that we only need to sum if we are running on multiple
  // processors, and we only need to take the square-root
  // if the value is nonzero.  There will in general be many
  // zeros for the inactive elements.

  // First sum the vector
  this->reduce_error(error_per_cell);

  // Compute the square-root of each component.
  START_LOG("std::sqrt()", "KellyErrorEstimator");
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = std::sqrt(error_per_cell[i]);
  STOP_LOG("std::sqrt()", "KellyErrorEstimator");
  
  //  STOP_LOG("flux_jumps()", "KellyErrorEstimator");
}







void
KellyErrorEstimator::attach_flux_bc_function (std::pair<bool,Real> fptr(const System& system,
									const Point& p,
									const std::string& var_name))
{
  assert (fptr != NULL);
  
  _bc_function = fptr;
}
