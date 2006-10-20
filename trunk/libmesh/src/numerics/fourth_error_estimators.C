// $Id: fourth_error_estimators.C,v 1.13 2006-10-20 20:31:22 roystgnr Exp $

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
#include "fourth_error_estimators.h"
#include "dof_map.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_clough.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
#include "system.h"

//-----------------------------------------------------------------
// ErrorEstimator implementations
#ifndef ENABLE_SECOND_DERIVATIVES

void LaplacianErrorEstimator::estimate_error (const System&,
					      ErrorVector&,
					      bool)
{
  std::cerr << "ERROR:  This functionalitry requires second-derivative support!"
	    << std::endl;
  error();
}

#else // defined (ENABLE_SECOND_DERIVATIVES)

void LaplacianErrorEstimator::estimate_error (const System& system,
					      ErrorVector& error_per_cell,
					      bool)
{
  START_LOG("laplacian_jump()", "LaplacianErrorEstimator");
  
  /*

  - e & f are global element ids
  
  Case (1.) Elements are at the same level, e<f
            Compute the laplacian jump on the face and
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
             Compute the laplacian jump on e's face and
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

  // Check for the use of component_mask
  this->convert_component_mask_to_scale();

  // Check for a valid component_scale
  if (!component_scale.empty())
    {
      if (component_scale.size() != n_vars)
	{
	  std::cerr << "ERROR: component_scale is the wrong size:"
		    << std::endl
		    << " component_scale.size()=" << component_scale.size()
		    << std::endl
		    << " n_vars=" << n_vars
		    << std::endl;
	  error();
	}
    }
  else
    {
      // No specified scaling.  Scale all variables by one.
      component_scale.resize (n_vars);
      std::fill (component_scale.begin(), component_scale.end(), 1.0);
    }
  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
	if (component_scale[var] == 0.0) continue;
      
      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);

      // Finite element objects for the same face from
      // different sides
      AutoPtr<FEBase> fe_e (FEBase::build (dim, fe_type));
      AutoPtr<FEBase> fe_f (FEBase::build (dim, fe_type));

      // Build an appropriate quadrature rule
      AutoPtr<QBase> qrule(fe_type.default_quadrature_rule(dim-1));

      // Tell the finite element for element e about the quadrature
      // rule.  The finite element for element f need not know about it
      fe_e->attach_quadrature_rule (qrule.get());
      
      // By convention we will always do the integration
      // on the face of element e.  Get its Jacobian values, etc..
      const std::vector<Real>&  JxW_face     = fe_e->get_JxW();
      const std::vector<Point>& qface_point  = fe_e->get_xyz();
//      const std::vector<Point>& face_normals = fe_e->get_normals();

      // The quadrature points on element f.  These will be computed
      // from the quadrature points on element e.
      std::vector<Point> qp_f;
      
      // The shape function second derivatives on elements e & f
      const std::vector<std::vector<RealTensor> > & d2phi_e =
		      fe_e->get_d2phi();
      const std::vector<std::vector<RealTensor> > & d2phi_f =
		      fe_f->get_d2phi();
      
      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices_e;
      std::vector<unsigned int> dof_indices_f;


      
      // Iterate over all the active elements in the mesh
      // that live on this processor.

      MeshBase::const_element_iterator       elem_it  =
		      mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end =
		      mesh.active_local_elements_end(); 

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* e = *elem_it;
	  const unsigned int e_id = e->id();
	  
	  // Loop over the neighbors of element e
	  for (unsigned int n_e=0; n_e<e->n_neighbors(); n_e++)
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
		    fe_e->reinit (e, n_e);

		    // Build the side
		    AutoPtr<Elem> side (e->build_side(n_e));

		    // Get the maximum h for this side
		    Real h_e, h_f;

                    if (dim == 1)
                      {
                        h_e = e->hmax();
                        h_f = f->hmax();
                      }
                    else
                      h_e = h_f = side->hmax();

		    // Get the DOF indices for the two elements
		    dof_map.dof_indices (e, dof_indices_e, var);
		    dof_map.dof_indices (f, dof_indices_f, var);

		    // The number of DOFS on each element
		    const unsigned int n_dofs_e = dof_indices_e.size();
		    const unsigned int n_dofs_f = dof_indices_f.size();

		    // The number of quadrature points
		    const unsigned int n_qp = qrule->n_points();

		    // Find the location of the quadrature points
		    // on element f
		    FEInterface::inverse_map (dim, fe_type, f,
					      qface_point, qp_f);

		    // Compute the shape functions on element f
		    // at the quadrature points of element e
		    fe_f->reinit (f, &qp_f);

		    // The error contribution from this face
		    Real error = 1.e-10;

		    
		    
		    // loop over the integration points on the face
		    for (unsigned int qp=0; qp<n_qp; qp++)
		      {
			// The solution laplacian from each element
			Number lap_e = 0., lap_f = 0.;
			
			// Compute the solution laplacian on element e
			for (unsigned int i=0; i<n_dofs_e; i++)
			  {
			  lap_e += d2phi_e[i][qp](0,0) *
				system.current_solution(dof_indices_e[i]);
                          if (dim > 1)
                            {
			  lap_e += d2phi_e[i][qp](1,1) *
				system.current_solution(dof_indices_e[i]);
                            }
                          if (dim > 2)
                            {
			  lap_e += d2phi_e[i][qp](2,2) *
				system.current_solution(dof_indices_e[i]);
                            }
			  }
			
			// Compute the solution gradient on element f
			for (unsigned int i=0; i<n_dofs_f; i++)
			  {
			  lap_f += d2phi_f[i][qp](0,0) *
				system.current_solution(dof_indices_f[i]);
                          if (dim > 1)
                            {
			  lap_f += d2phi_f[i][qp](1,1) *
				system.current_solution(dof_indices_f[i]);
                            }
                          if (dim > 2)
                            {
			  lap_f += d2phi_f[i][qp](2,2) *
				system.current_solution(dof_indices_f[i]);
                            }
			  }
			

			// The flux jump at the face 
			const Number jump = lap_e - lap_f;

			// The flux jump squared
#ifndef USE_COMPLEX_NUMBERS
			const Real jump2 = jump*jump;
#else
			const Real jump2 = std::norm(jump);
#endif

			// Integrate the error on the face
			error += JxW_face[qp]*jump2;
			
		      } // End quadrature point loop

		    // Add the error contribution to elements e & f
		    error_per_cell[e_id] += error*h_e*component_scale[var];
		    error_per_cell[f_id] += error*h_f*component_scale[var];
		    
		    
		  }
	      } // if (e->neigbor(n_e) != NULL)
	  
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
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = sqrt(error_per_cell[i]);
  
  STOP_LOG("laplacian_jump()", "LaplacianErrorEstimator");
}

#endif // defined (ENABLE_SECOND_DERIVATIVES)
