// $Id: jump_error_estimator.C,v 1.2 2006-10-26 17:15:59 roystgnr Exp $

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
#include "jump_error_estimator.h"
#include "dof_map.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
#include "system.h"

#include "dense_vector.h"

//-----------------------------------------------------------------
// JumpErrorEstimator implementations
void JumpErrorEstimator::initialize (const System&,
				     ErrorVector&,
				     bool)
{
}



void JumpErrorEstimator::estimate_error (const System& system,
					 ErrorVector& error_per_cell,
					 bool estimate_parent_error)
{
  START_LOG("estimate_error()", "JumpErrorEstimator");
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

  // Declare a vector of floats which is as long as
  // error_per_cell above, and fill with zeros.  This vector will be
  // used to keep track of the number of edges (faces) on each active
  // element which are either:
  // 1) an internal edge
  // 2) an edge on a Neumann boundary for which a boundary condition
  //    function has been specified.
  // The error estimator can be scaled by the number of flux edges (faces)
  // which the element actually has to obtain a more uniform measure
  // of the error.  Use floats instead of ints since in case 2 (above)
  // f gets 1/2 of a flux face contribution from each of his
  // neighbors
  std::vector<float> n_flux_faces (error_per_cell.size());
  
  // Check for the use of component_mask
  this->convert_component_mask_to_scale();

  // Check for a valid component_scale
  if (!component_scale.empty())
    {
      if (component_scale.size() != n_vars)
	{
	  std::cerr << "ERROR: component_scale is the wrong size:"
		    << std::endl
		    << " component_scale.scale()=" << component_scale.size()
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
  for (var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
	if (component_scale[var] == 0.0) continue;

      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);
      
      // Finite element objects for the same face from
      // different sides
      fe_fine = FEBase::build (dim, fe_type);
      fe_coarse = FEBase::build (dim, fe_type);

      // Build an appropriate Gaussian quadrature rule
      QGauss qrule (dim-1, fe_type.default_quadrature_order());

      // Tell the finite element for the fine element about the quadrature
      // rule.  The finite element for the coarse element need not know about it
      fe_fine->attach_quadrature_rule (&qrule);
      
      // By convention we will always do the integration
      // on the face of element e.  We'll need its Jacobian values and
      // physical point locations, at least
      fe_fine->get_JxW();
      fe_fine->get_xyz();

      // Our derived classes may want to do some initialization here
      this->initialize(system, error_per_cell, estimate_parent_error);
      
      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices_fine;
      std::vector<unsigned int> dof_indices_coarse;


      
      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* e = *elem_it;
	  const unsigned int e_id = e->id();

          // See if the parent of element e has been examined yet;
          // if not, we may want to compute the estimator on it
          const Elem* parent = e->parent();

          // We only can compute and only need to compute on
          // parents with all active children
          bool compute_on_parent = true;
          if (!parent || !estimate_parent_error)
            compute_on_parent = false;
          else
            for (unsigned int c=0; c != parent->n_children(); ++c)
              if (!parent->child(c)->active())
                compute_on_parent = false;
             
          if (compute_on_parent &&
              !error_per_cell[parent->id()])
	    {
              // Compute a projection onto the parent
              DenseVector<Number> Uparent;
              FEBase::coarsened_dof_values(*(system.solution),
                                           dof_map, parent, Uparent,
                                           var, false);

	      // Loop over the neighbors of the parent
	      for (unsigned int n_p=0; n_p<parent->n_neighbors(); n_p++)
                {
	          if (parent->neighbor(n_p) != NULL) // parent has a neighbor here
		    {
                      // Find the active neighbors in this direction
                      std::vector<const Elem*> active_neighbors;
                      parent->active_family_tree_by_neighbor(active_neighbors,
                                                             parent);
                      // Compute the flux to each active neighbor
                      for (unsigned int a=0; 
                           a != active_neighbors.size(); ++a)
                        {
                          const Elem *f = active_neighbors[a];
                          if (f->level() >= parent->level())
                            {
                              fine_elem = f;
                              coarse_elem = parent;
                              Ucoarse = Uparent;

		              dof_map.dof_indices (fine_elem, dof_indices_fine, var);
		              const unsigned int n_dofs_fine = dof_indices_fine.size();
                              Ufine.resize(n_dofs_fine);
			
			      for (unsigned int i=0; i<n_dofs_fine; i++)
			        Ufine(i) = system.current_solution(dof_indices_fine[i]);
                              this->reinit_sides();
                              this->internal_side_integration();

                              error_per_cell[fine_elem->id()] += fine_error;
                              error_per_cell[coarse_elem->id()] += coarse_error;

                              // Keep track of the number of internal flux
                              // sides found on each element
                              n_flux_faces[fine_elem->id()]++;
                              n_flux_faces[coarse_elem->id()] += this->coarse_n_flux_faces_increment();
                            }
                        }
		    }
		  else if (integrate_boundary_sides)
		    {
                      fine_elem = parent;
                      Ufine = Uparent;

                      // Reinitialize shape functions on the fine element side
                      fe_fine->reinit (fine_elem, fine_side);

                      if (this->boundary_side_integration())
                        {
                          error_per_cell[fine_elem->id()] += fine_error;
                          n_flux_faces[fine_elem->id()]++;
                        }
                    } 
		}
	    }

          // If we do any more flux integration, e will be the fine element
          fine_elem = e;

	  // Loop over the neighbors of element e
	  for (unsigned int n_e=0; n_e<e->n_neighbors(); n_e++)
	    {
              fine_side = n_e;

	      if (e->neighbor(n_e) != NULL) // e is not on the boundary
		{
		  const Elem* f           = e->neighbor(n_e);
		  const unsigned int f_id = f->id();

		  // Compute flux jumps if we are in case 1 or case 2.
		  if ((f->active() && (f->level() == e->level()) && (e_id < f_id))
		      || (f->level() < e->level()))
		    {		    
                      // f is now the coarse element
                      coarse_elem = f;

		      // Get the DOF indices for the two elements
		      dof_map.dof_indices (fine_elem, dof_indices_fine, var);
		      dof_map.dof_indices (coarse_elem, dof_indices_coarse, var);

		      // The number of DOFS on each element
		      const unsigned int n_dofs_fine = dof_indices_fine.size();
		      const unsigned int n_dofs_coarse = dof_indices_coarse.size();
                      Ufine.resize(n_dofs_fine);
                      Ucoarse.resize(n_dofs_coarse);

		      // The local solutions on each element
		      for (unsigned int i=0; i<n_dofs_fine; i++)
			Ufine(i) = system.current_solution(dof_indices_fine[i]);
		      for (unsigned int i=0; i<n_dofs_coarse; i++)
			Ucoarse(i) = system.current_solution(dof_indices_coarse[i]);
			
                      this->reinit_sides();
                      this->internal_side_integration();

                      error_per_cell[fine_elem->id()] += fine_error;
                      error_per_cell[coarse_elem->id()] += coarse_error;

                      // Keep track of the number of internal flux
                      // sides found on each element
                      n_flux_faces[fine_elem->id()]++;
                      n_flux_faces[coarse_elem->id()] += this->coarse_n_flux_faces_increment();
		    } // end if (case1 || case2)
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
		  if (integrate_boundary_sides)
		    {
                      // Reinitialize shape functions on the fine element side
                      fe_fine->reinit (fine_elem, fine_side);

		      // Get the DOF indices 
		      dof_map.dof_indices (fine_elem, dof_indices_fine, var);

		      // The number of DOFS on each element
		      const unsigned int n_dofs_fine = dof_indices_fine.size();
                      Ufine.resize(n_dofs_fine);

                      for (unsigned int i=0; i<n_dofs_fine; i++)
                        Ufine(i) = system.current_solution(dof_indices_fine[i]);

                      if (this->boundary_side_integration())
                        {
                          error_per_cell[fine_elem->id()] += fine_error;
                          n_flux_faces[fine_elem->id()]++;
                        }
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

  // First sum the vector of estimated error values
  this->reduce_error(error_per_cell);

  // Compute the square-root of each component.
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = std::sqrt(error_per_cell[i]);


  if (this->scale_by_n_flux_faces)
    {
      // Sum the vector of flux face counts
      this->reduce_error(n_flux_faces);

      // Sanity check: Make sure the number of flux faces is
      // always an integer value
#ifdef DEBUG
      for (unsigned int i=0; i<n_flux_faces.size(); ++i)
	assert (n_flux_faces[i] == static_cast<float>(static_cast<unsigned int>(n_flux_faces[i])) );
#endif
  
      // Scale the error by the number of flux faces for each element
      for (unsigned int i=0; i<n_flux_faces.size(); ++i)
	{
	  if (n_flux_faces[i] == 0.0) // inactive or non-local element
	    continue;
      
	  //std::cout << "Element " << i << " has " << n_flux_faces[i] << " flux faces." << std::endl;
	  error_per_cell[i] /= static_cast<Real>(n_flux_faces[i]); 
	}
    }
  
  STOP_LOG("estimate_error()", "JumpErrorEstimator");
}



void
JumpErrorEstimator::reinit_sides ()
{
  // The master quadrature point locations on the coarse element
  std::vector<Point> qp_coarse;

  // Reinitialize shape functions on the fine element side
  fe_fine->reinit (fine_elem, fine_side);

  // Get the physical locations of the fine element quadrature points
  std::vector<Point> qface_point = fe_fine->get_xyz();

  // Find their locations on the coarse element
  FEInterface::inverse_map (coarse_elem->dim(), fe_coarse->get_fe_type(),
                            coarse_elem, qface_point, qp_coarse);

  // Calculate the coarse element shape functions at those locations
  fe_coarse->reinit (coarse_elem, &qp_coarse);
}



float JumpErrorEstimator::coarse_n_flux_faces_increment ()
{
  // Keep track of the number of internal flux sides found on each
  // element
  unsigned int dim = coarse_elem->dim();

  const unsigned int divisor =
    1 << (dim-1)*(fine_elem->level() - coarse_elem->level());

  // With a difference of n levels between fine and coarse elements,
  // we compute a fractional flux face for the coarse element by adding:
  // 1/2^n in 2D
  // 1/4^n in 3D
  // each time.  This code will get hit 2^n times in 2D and 4^n
  // times in 3D so that the final flux face count for the coarse
  // element will be an integer value.

  return 1.0 / static_cast<Real>(divisor);
}
