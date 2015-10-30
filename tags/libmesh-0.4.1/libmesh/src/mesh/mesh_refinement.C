// $Id: mesh_refinement.C,v 1.29 2003-11-05 22:26:45 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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

#include "mesh_base.h"
#include "mesh_refinement.h"
#include "elem.h"
#include "libmesh_logging.h"



//-----------------------------------------------------------------
// Mesh refinement methods
MeshRefinement::MeshRefinement (MeshBase& m) :
  _mesh(m)
{
}



MeshRefinement::~MeshRefinement ()
{
  this->clear();  
}



void MeshRefinement::clear ()
{
  _new_nodes_map.clear();
  _unused_elements.clear();
}



Node* MeshRefinement::add_point (const Point& p, const unsigned int key)
{
  // Look for the key in the multimap  
  std::pair<map_type::iterator, map_type::iterator>
    pos = _new_nodes_map.equal_range(key);
  
      
  while (pos.first != pos.second) 
    if (p == *pos.first->second)
      return pos.first->second;      
    else      
      ++pos.first;
    

  // If we get here pos.first == pos.second.
  assert (pos.first == pos.second); // still not found
                                    // so we better add it
  Node* node = _mesh.add_point (p);

  assert (node != NULL);

  // Add the node to the map.  In the case of the
  // std::multimap use pos.first as a hint for where to put it
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
  _new_nodes_map.insert(std::make_pair(key, node));
#else
  _new_nodes_map.insert(pos.first, std::make_pair(key, node));
#endif			    

  // Set the key for this node
  node->set_key() = key;
  
  // Return the address of the new node
  return node;
}



Elem* MeshRefinement::add_elem (Elem* elem)
{
  assert (elem != NULL);
  
  // If the unused_elements has any iterators from
  // old elements, take the first one
  if (!_unused_elements.empty())
    {
      std::vector<Elem*>::iterator it = _unused_elements.front();

      *it = elem;

      _unused_elements.pop_front();
    }

  // Otherwise, use the conventional add method
  else
    {
      _mesh.add_elem (elem);
    }

  return elem;
}



void MeshRefinement::update_nodes_map ()
{
  START_LOG("update_nodes_map()", "MeshRefinement");

  _new_nodes_map.clear();
  
  std::vector<Node*> & nodes = _mesh.get_nodes();

  std::vector<Node*>::iterator       it  = nodes.begin();
  const std::vector<Node*>::iterator end = nodes.end();

  for (; it != end; ++it)
    {
      Node* node = *it;

      assert (node != NULL);
      
      // Add the node to the map.
      _new_nodes_map.insert(std::make_pair(node->key(), node));
    }

  STOP_LOG("update_nodes_map()", "MeshRefinement");
} 



void MeshRefinement::refine_and_coarsen_elements (const bool maintain_level_one)
{
  assert (_mesh.mesh_dimension() != 1);

  bool satisfied = false;


  // Possibly clean up the refinement flags from
  // a previous step
  elem_iterator       elem_it (_mesh.elements_begin());
  const elem_iterator elem_end(_mesh.elements_end());

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Set refinement flag to DO_NOTHING if the
      // element isn't active
      if ( !(*elem_it)->active())
	(*elem_it)->set_refinement_flag(Elem::DO_NOTHING);

      // This might be left over from the last step
      if ((*elem_it)->refinement_flag() == Elem::JUST_REFINED)
	(*elem_it)->set_refinement_flag(Elem::DO_NOTHING);
    }

  
  /**
   * Repeat until coarsening & refinement flags jive
   */
  do
    {
      const bool coarsening_satisfied =
	this->make_coarsening_compatible(maintain_level_one);
      
      const bool refinement_satisfied =
	this->make_refinement_compatible(maintain_level_one);

      const bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();// &&
// 	!this->limit_level_mismatch_at_node(1);
      
      satisfied = (coarsening_satisfied &&
		   refinement_satisfied &&
		   smoothing_satisfied);
    }
  while (!satisfied);

  /**
   * First coarsen the flagged elements.  This
   * will free space.
   */
  this->coarsen_elements ();

  /**
   * Now refine the flagged elements.  This will
   * take up some space, maybe more than what was freed.
   */
  this->refine_elements();
  
  /**
   * Finally, the new mesh needs to be prepared for use
   */
  _mesh.prepare_for_use ();
}



bool MeshRefinement::make_coarsening_compatible(const bool maintain_level_one)
{
  START_LOG ("make_coarsening_compatible()", "MeshRefinement");
  
  /**
   * Unless we encounter a specific situation level-one
   * will be satisfied after executing this loop just once
   */
  bool level_one_satisfied = true;

  /**
   * Unless we encounter a specific situation we will be compatible
   * with any selected refinement flags
   */
  bool compatible_with_refinement = true;

  /**
   * find the maximum level in the mesh
   */
  unsigned int max_level = 0;
    
  /**
   * First we look at all the active level-0 elements.  Since it doesn't make
   * sense to coarsen them we must un-set their coarsen flags if
   * they are set.
   */
  for (unsigned int e=0; e<_mesh.n_elem(); e++)
    {      
      assert (_mesh.elem(e) != NULL);
      max_level = std::max(max_level, _mesh.elem(e)->level());
      
      if (_mesh.elem(e)->active()      &&
	  (_mesh.elem(e)->level() == 0) &&
	  (_mesh.elem(e)->refinement_flag() == Elem::COARSEN))
	_mesh.elem(e)->set_refinement_flag(Elem::DO_NOTHING);
    }

  /**
   * if there are no refined elements then
   * there is no work for us to do
   */
  if (max_level == 0)
    {
      STOP_LOG ("make_coarsening_compatible()", "MeshRefinement");
      return compatible_with_refinement;
    }

  
  /**
   * Loop over all the active elements.  If an element is marked
   * for coarsening we better check its neighbors.  If ANY of these neighbors
   * are marked for refinement AND are at the same level then there is a
   * conflict.  By convention refinement wins, so we un-mark the element for
   * coarsening.  Level-one would be violated in this case so we need to re-run
   * the loop.
   */
  if (maintain_level_one)
    {
         
    repeat:
      level_one_satisfied = true;
  
      do
	{
	  level_one_satisfied = true;
	  
	  active_elem_iterator       el    (_mesh.elements_begin());
	  const active_elem_iterator end_el(_mesh.elements_end());
	  
	  for (; el != end_el; ++el)
	    {
	      Elem* elem = *el;
	      
	      if (elem->refinement_flag() == Elem::COARSEN) // If the element is active and 
		// the coarsen flag is set
		{
		  const unsigned int my_level = elem->level();
		  
		  for (unsigned int n=0; n<elem->n_neighbors(); n++)
		    if (elem->neighbor(n) != NULL)     // I have a neighbor
		      if (elem->neighbor(n)->active()) // and it is active
			{
			  const Elem* neighbor = elem->neighbor(n);
			  
			  if ((neighbor->level() == my_level) &&
			      (neighbor->refinement_flag() == Elem::REFINE)) // the neighbor is at my level
	        	                                                     // and wants to be refined
			    {
			      elem->set_refinement_flag(Elem::DO_NOTHING);
			      level_one_satisfied = false;
			    }
			}
		      else // I have a neighbor and it is not active. That means it has children.
			{  // While it _may_ be possible to coarsen us if all the children of
			   // that element want to be refined, it is impossible to know at this
			   // stage.  Forget about it for the moment...  This can be handled in
			   // two steps.
			  elem->set_refinement_flag(Elem::DO_NOTHING);
			  level_one_satisfied = false;
			}
		}
	    }      
	}
      while (!level_one_satisfied);
      
    } // end if (maintian_level_one)
  
  /**
   * Next we look at all of the inactive cells.
   * If there is an inactive cell with all of its children
   * wanting to be unrefined then the element is a candidate
   * for unrefinement.  If all the children don't
   * all want to be unrefined then ALL of them need to have their
   * unrefinement flags cleared. 
   */
  for (int level=(max_level); level >= 0; level--)
    for (unsigned int e=0; e<_mesh.n_elem(); e++)
      if (!_mesh.elem(e)->active() &&
	  (_mesh.elem(e)->level() == static_cast<unsigned int>(level)))
	{
	  /**
	   * right now the element hasn't been disqualified
	   * as a candidate for unrefinement
	   */
	  bool is_a_candidate = true;
	      
	  for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
	    if ((_mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN) ||
		!_mesh.elem(e)->child(c)->active() )
	      is_a_candidate = false;
	      
	  if (!is_a_candidate)
	    {
	      _mesh.elem(e)->set_refinement_flag(Elem::DO_NOTHING);
		  
	      for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
		if (_mesh.elem(e)->child(c)->refinement_flag() == Elem::COARSEN)
		  {
		    level_one_satisfied = false;
		    _mesh.elem(e)->child(c)->set_refinement_flag(Elem::DO_NOTHING);
		  }
	    }
	}
      
  if (!level_one_satisfied && maintain_level_one) goto repeat;
  
  
  /**
   * If all the children of a parent are set to be coarsened
   * then flag the parent so that they can kill thier kids...
   */
  for (unsigned int e=0; e<_mesh.n_elem(); e++)
    if (!_mesh.elem(e)->active())
      {
	/**
	 * Presume all the children are flagged for coarsening
	 * and then look for a contradiction
	 */
	bool all_children_flagged_for_coarsening = true;
	
	for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
	  if (_mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN)
	    all_children_flagged_for_coarsening = false;
	
	if (all_children_flagged_for_coarsening)
	  _mesh.elem(e)->set_refinement_flag(Elem::COARSEN);
      }
	
  STOP_LOG ("make_coarsening_compatible()", "MeshRefinement");
  
  return compatible_with_refinement;
}



bool MeshRefinement::make_refinement_compatible(const bool maintain_level_one)
{
  START_LOG ("make_refinement_compatible()", "MeshRefinement");
  
  // Unless we encounter a specific situation level-one
  // will be satisfied after executing this loop just once
  bool level_one_satisfied = true;

  // Unless we encounter a specific situation we will be compatible
  // with any selected coarsening flags
  bool compatible_with_coarsening = true;

#ifdef DEBUG
  
  // It only makes sense to refine active elements
  for (unsigned int e=0; e<_mesh.n_elem(); e++)
    if (_mesh.elem(e)->refinement_flag() == Elem::REFINE)
      assert (_mesh.elem(e)->active());
  
#endif
  

  /**
   * This loop enforces the level-1 rule.  We should only
   * execute it if the user indeed wants level-1 satisfied!
   */
  if (maintain_level_one)
    {  
      do
	{
	  level_one_satisfied = true;
	  
	  for (unsigned int e=0; e<_mesh.n_elem(); e++)
	    if (_mesh.elem(e)->active() &&                
		(_mesh.elem(e)->refinement_flag() == Elem::REFINE))  // If the element is active and the
	                                                            // refinement flag is set    
	      {
		const unsigned int my_level = _mesh.elem(e)->level();
	    
		for (unsigned int side=0; side<_mesh.elem(e)->n_sides(); side++)
		  if (_mesh.elem(e)->neighbor(side) != NULL)     // I have a neighbor
		    if (_mesh.elem(e)->neighbor(side)->active()) // and it is active
		      {
			Elem* neighbor = _mesh.elem(e)->neighbor(side);
			
			/**
			 * Case 1:  The neighbor is at the same level I am.
			 *        1a: The neighbor will be refined       -> NO PROBLEM
			 *        1b: The neighbor won't be refined      -> NO PROBLEM
			 *        1c: The neighbor wants to be coarsened -> PROBLEM
			 */
			if (neighbor->level() == my_level)
			  {
			    if (neighbor->refinement_flag() == Elem::COARSEN)
			      {
				neighbor->set_refinement_flag(Elem::DO_NOTHING);
				compatible_with_coarsening = false;
				level_one_satisfied = false;
			      }
			  }
			
			/**
			 * Case 2: The neighbor is one level lower than I am.
			 *         The neighbor thus MUST be refined to satisfy
			 *         the level-one rule, regardless of whether it
			 *         was originally flagged for refinement. If it
			 *         wasn't flagged already we need to repeat
			 *         this process.
			 */
			else if ((neighbor->level()+1) == my_level)
			  {
			    if (neighbor->refinement_flag() != Elem::REFINE)
			      {
				neighbor->set_refinement_flag(Elem::REFINE);
				level_one_satisfied = false; 
			      }
			  }
#ifdef DEBUG		
			/**
			 * Sanity check. We should never get into a
			 * case when our neighbot is more than one
			 * level away.
			 */
			else if ((neighbor->level()+1) < my_level)
			  {
			    error();
			  }
			
			/**
			 * Note that the only other possibility is that the
			 * neighbor is already refined, in which case it isn't
			 * active and we should never get here.
			 */
			else
			  {
			    error();
			  }
#endif		
		      } 
	      }
	}
      
      while (!level_one_satisfied);
    } // end if (maintain_level_one)

  
  STOP_LOG ("make_refinement_compatible()", "MeshRefinement");
  
  return compatible_with_coarsening;
}



void MeshRefinement::coarsen_elements ()
{
  START_LOG ("coarsen_elements()", "MeshRefinement");

  // Clear the unused_elements data structure.
  // The elements have been packed since it was built,
  // so there are _no_ unused elements.  We cannot trust
  // any iterators currently in this data structure.
  _unused_elements.clear();

  std::vector<Elem*>& elements = _mesh.get_elem();

  // Loop over the elements.   
  for (std::vector<Elem*>::iterator it = elements.begin();
       it != elements.end(); ++it)
    {
      Elem* elem = *it;

      assert (elem != NULL);
      
      if (elem->refinement_flag() == Elem::COARSEN)
	// active elements flagged for coarsening will
	// truly be deleted
	if (elem->active()) 
	  {
	    // Huh?  no level-0 element should be active
	    // and flagged for coarsening.
	    assert (elem->level() != 0);

	    // Remove this element from any neighbor
	    // lists that point to it.
	    elem->nullify_neighbors();

	    // Remove any boundary information associated
	    // with this element
	    _mesh.boundary_info.remove (elem);
	    
	    // Delete the element, make sure the corresponding
	    // entry in the vector is set to NULL
	    delete elem;

	    *it = NULL;
	    
	    // Add this iterator to the _unused_elements
	    // data structure so we might fill it.
	    _unused_elements.push_back (it);
	  }
        // inactive elements flagged for coarsening
        // will become active
	else
	  {
	    elem->coarsen();
	    assert (elem->active());
	  }
    }
  
  STOP_LOG ("coarsen_elements()", "MeshRefinement");  
}



void MeshRefinement::refine_elements ()
{
  // Update the _new_nodes_map so that elements can
  // find nodes to connect to.
  this->update_nodes_map ();

  START_LOG ("refine_elements()", "MeshRefinement");
  
  // Careful!  We can't use iterators here.  When
  // an element is refined reallocation in the elements
  // vector may be forced, which might invalidate
  // any iterators we have.
  const unsigned int orig_n_elem = _mesh.n_elem();

  std::vector<Elem*>& elements = _mesh.get_elem();
  
  // Loop over the elements in the mesh, refine
  // those whose flags are set to Elem::REFINE
  for (unsigned int e=0; e<orig_n_elem; e++)
    if (elements[e] != NULL)
      if (elements[e]->refinement_flag() == Elem::REFINE)
	elements[e]->refine (*this);

  this->clear();
  
  STOP_LOG ("refine_elements()", "MeshRefinement");
}



void MeshRefinement::uniformly_refine (unsigned int n)
{
  // Refine n times
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // Clean up the refinement flags
      this->clean_refinement_flags();
      
      // Flag all the active elements for refinement.       
      active_elem_iterator       elem_it (_mesh.elements_begin());
      const active_elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_refinement_flag(Elem::REFINE);

      // Refine all the elements we just flagged.
      this->refine_elements();
    }
  
  // Finally, the new mesh needs to be prepared for use
  _mesh.prepare_for_use ();
}


#endif
