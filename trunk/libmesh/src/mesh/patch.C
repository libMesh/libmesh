
// $Id: patch.C,v 1.1 2007-01-18 22:24:47 roystgnr Exp $

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
#include "libmesh_logging.h"
#include "patch.h"
#include "elem.h"



//-----------------------------------------------------------------
// Patch implementations
void Patch::add_face_neighbors()
{
  std::set<const Elem *> new_neighbors;
    
  // Loop over all the elements in the patch
  std::set<const Elem*>::const_iterator       it  = this->begin();
  const std::set<const Elem*>::const_iterator end = this->end();

  for (; it != end; ++it)
    {
      const Elem* elem = *it;
      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor(s) != NULL)        // we have a neighbor on this side
          {
	    const Elem* neighbor = elem->neighbor(s);
	
	    if (neighbor->active())           // ... and that neighbor is active
	      new_neighbors.insert (neighbor);// ... then add it to the patch
	  
	    else                              // ... the neighbor is *not* active,
	      {                               // ... so add *all* neighboring
                                              // active children to the patch
	        std::vector<const Elem*> active_neighbor_children;

	        neighbor->active_family_tree_by_neighbor
                  (active_neighbor_children, elem);

		new_neighbors.insert(active_neighbor_children.begin(),
                                     active_neighbor_children.end());
	      }
          }
    }

  this->insert(new_neighbors.begin(), new_neighbors.end());
}



void Patch::add_local_face_neighbors()
{
  std::set<const Elem *> new_neighbors;
    
  // Loop over all the elements in the patch
  std::set<const Elem*>::const_iterator       it  = this->begin();
  const std::set<const Elem*>::const_iterator end = this->end();

  for (; it != end; ++it)
    {
      const Elem* elem = *it;
      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor(s) != NULL)        // we have a neighbor on this side
          {
	    const Elem* neighbor = elem->neighbor(s);
	
	    if (neighbor->active())           // ... and that neighbor is active
	      if (neighbor->processor_id() ==
	          libMesh::processor_id())    // ... and belongs to this processor
	        new_neighbors.insert (neighbor);// ... then add it to the patch
	  
	    else                              // ... the neighbor is *not* active,
	      {                               // ... so add *all* neighboring
                                              // active, local children to the patch
	        std::vector<const Elem*> active_neighbor_children;

	        neighbor->active_family_tree_by_neighbor
                  (active_neighbor_children, elem);

	        for (unsigned int c=0; c<active_neighbor_children.size(); c++)
	          if (active_neighbor_children[c]->processor_id() == libMesh::processor_id())
		    new_neighbors.insert(active_neighbor_children[c]);
	      }
          }
    }

  this->insert(new_neighbors.begin(), new_neighbors.end());
}
    


void Patch::add_vertex_neighbors()
{
  std::cerr << "ERROR: Patch::add_vertex_neighbors() is not yet implemented"
	    << std::endl;
  error();
}


  
void Patch::add_local_vertex_neighbors()
{
  std::cerr << "ERROR: Patch::add_local_vertex_neighbors() is not yet implemented"
	    << std::endl;
  error();
}
  


void Patch::build_around_element (const Elem* e0,
                                  const unsigned int target_patch_size,
                                  PMF patchtype)
{
  START_LOG("build_around_element()", "Patch");
  
  // Make sure we are building a patch for an active element.
  assert (e0 != NULL);
  assert (e0->active());
  // Make sure we are either starting with a local element or
  // requesting a nonlocal patch
  assert ((patchtype != &Patch::add_local_face_neighbors &&
           patchtype != &Patch::add_local_vertex_neighbors) ||
           e0->processor_id() == libMesh::processor_id());
  
  // First clear the current set, then add the element of interest.
  this->clear();
  this->insert (e0);

  // Repeatedly add the neighbors of the elements in the patch until
  // the target patch size is met
  while (this->size() < target_patch_size)
    {
      // It is possible that the target patch size is larger than the number
      // of elements that can be added to the patch.  Since we don't
      // have access to the Mesh object here, the only way we can
      // detect this case is by detecting a "stagnant patch," i.e. a
      // patch whose size does not increase after adding face neighbors
      const unsigned int old_patch_size = this->size();
      
      (this->*patchtype)();
      
      // Check for a "stagnant" patch
      if (this->size() == old_patch_size)
	{
	  std::cerr << "WARNING: stagnant patch of "
		    << this->size() << " elements."
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
    std::set<const Elem*>::const_iterator       it  = this->begin();
    const std::set<const Elem*>::const_iterator end = this->end();
    
    for (; it != end; ++it)
      {
	// Convenience.  Keep the syntax simple.
	const Elem* elem = *it;

	assert (elem->active());
        if ((patchtype == &Patch::add_local_face_neighbors ||
             patchtype == &Patch::add_local_vertex_neighbors))
	  assert (elem->processor_id() == libMesh::processor_id());
      }
  }
#endif

  STOP_LOG("build_around_element()", "Patch");
}
