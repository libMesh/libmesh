// $Id: mesh_refinement_smoothing.C,v 1.8 2004-11-08 00:11:05 jwpeterson Exp $

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



// Local includes
#include "libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef ENABLE_AMR

#include "mesh_refinement.h"
#include "mesh_base.h"
#include "elem.h"



//-----------------------------------------------------------------
// Mesh refinement methods
bool MeshRefinement::limit_level_mismatch_at_node (const unsigned int max_mismatch)
{
  bool flags_changed = false;


  // Vector holding the maximum element level that touches a node.
  std::vector<unsigned char> max_level_at_node (_mesh.n_nodes(), 0);


  // Loop over all the active elements & fill the vector
  {
//     const_active_elem_iterator       elem_it (_mesh.const_elements_begin());
//     const const_active_elem_iterator elem_end(_mesh.const_elements_end());

    MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 
    
    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;	
	const unsigned char elem_level =
	  elem->level() + ((elem->refinement_flag() == Elem::REFINE) ? 1 : 0);

	// Set the max_level at each node
	for (unsigned int n=0; n<elem->n_nodes(); n++)
	  {
	    const unsigned int node_number = elem->node(n);

	    assert (node_number < max_level_at_node.size());
	    
	    max_level_at_node[node_number] =
	      std::max (max_level_at_node[node_number], elem_level);
	  }
      }     
  }


  // Now loop over the active elements and flag the elements
  // who violate the requested level mismatch
  {
//     active_elem_iterator       elem_it (_mesh.elements_begin());
//     const active_elem_iterator elem_end(_mesh.elements_end());

    MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 
    
    for (; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;	
	const unsigned int elem_level = elem->level();
	
	// Skip the element if it is already flagged
	if (elem->refinement_flag() == Elem::REFINE)
	  continue;

	// Loop over the nodes, check for possible mismatch
	for (unsigned int n=0; n<elem->n_nodes(); n++)
	  {
	    const unsigned int node_number = elem->node(n);

	    // Flag the element for refinement if it violates
	    // the requested level mismatch
	    if ( (elem_level + max_mismatch) < max_level_at_node[node_number])
	      {
		elem->set_refinement_flag (Elem::REFINE);
		flags_changed = true;
	      }
	  }
      }     
  }
  
  return flags_changed;
}




bool MeshRefinement::eliminate_unrefined_patches ()
{
  bool flags_changed = false;
  
  
//   active_elem_iterator       elem_it (_mesh.elements_begin());
//   const active_elem_iterator elem_end(_mesh.elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;
      const unsigned int my_level = elem->level();
      bool flag_me = true;


      // Skip the element if it is already flagged
      if (elem->refinement_flag() == Elem::REFINE)
	continue;
      
      // Check all the element neighbors
      for (unsigned int n=0; n<elem->n_neighbors(); n++)
	{
	  // Quit if the element is on the boundary
	  if (elem->neighbor(n) == NULL)
	    {
	      flag_me = false;
	      break;
	    }

	  // Can't happen, would have broken the loop!
	  assert (elem->neighbor(n) != NULL);

	  // If the neighbor is active and will not be refined
	  // we do not need to refine ourselves.
	  if ((elem->neighbor(n)->level() == my_level) &&
	      (elem->neighbor(n)->active()) &&
	      (elem->neighbor(n)->refinement_flag() != Elem::REFINE))
	    flag_me = false;

	  // Otherwise the neighbor is inactive (i.e. already refined)
	  // or is active and will be refined.  In either case this
	  // neighb
	}

      if (flag_me)
	{
	  elem->set_refinement_flag(Elem::REFINE);
	  flags_changed = true;
	}      
    }

  return flags_changed;
}


#endif
