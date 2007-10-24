// $Id: mesh_refinement.C,v 1.12 2003-05-01 18:44:31 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef ENABLE_AMR

#include "mesh_refinement.h"
#include "mesh_base.h"
#include "elem.h"



//-----------------------------------------------------------------
// Mesh refinement methods
MeshRefinement::MeshRefinement (MeshBase& m) :
  mesh(m)
{
}



MeshRefinement::~MeshRefinement ()
{
  clear();  
}



void MeshRefinement::clear ()
{
  new_nodes.clear();
  unused_elements.clear();
}



Node* MeshRefinement::add_point(const Point& p)
{
  // if the new_nodes data structure is empty we need to fill it.
  if (new_nodes.empty())
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      new_nodes.insert(std::pair<unsigned int,
		                 unsigned int>(mesh.point(n).key(),n));

  
  unsigned int n;
  
  std::pair<std::multimap<unsigned int, unsigned int>::iterator,
            std::multimap<unsigned int, unsigned int>::iterator>
    pos = new_nodes.equal_range(p.key());
  
      
  while (pos.first != pos.second) 
    {                              
      if (mesh.point(pos.first->second) == p)
	{
	  n = pos.first->second;
	  
	  break;
	}
      ++pos.first;
    }
  
  if (pos.first == pos.second) // still not found
    {                          // so we better add it
      n = mesh.n_nodes();
      
      new_nodes.insert(pos.first,
		       std::pair<unsigned int,
		                 unsigned int>(p.key(),n));
    }
    
  
  // Return the address of the new node
  return mesh.add_point(p,n);
}



unsigned int MeshRefinement::new_element_number()
{
  unsigned int n=0;
  
  if (unused_elements.empty())      
    n = mesh.n_elem();
  else
    {
      n = *unused_elements.begin();
      unused_elements.erase(unused_elements.begin());
    }

  return n;
}



void MeshRefinement::update_unused_elements()
{
  /**
   * Update the unused_elements data structure
   */
  
  // clear and start from scratch
  unused_elements.clear();
  
  // find NULL elements and add the index to the database
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e) == NULL)
      unused_elements.insert(e);
}



void MeshRefinement::refine_and_coarsen_elements (const bool maintain_level_one)
{
  assert (mesh.mesh_dimension() != 1);

  bool satisfied = false;


  elem_iterator       elem_it (mesh.elements_begin());
  const elem_iterator elem_end(mesh.elements_end());

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
      const bool coarsening_satisfied = make_coarsening_compatible(maintain_level_one);
      const bool refinement_satisfied = make_refinement_compatible(maintain_level_one);
      
      if (coarsening_satisfied && refinement_satisfied)
	satisfied = true;
    }
  while (!satisfied);

  /**
   * First coarsen the flagged elements.  This
   * will free space.
   */
  coarsen_elements ();

  /**
   * Find any NULL elements from the coarsening step.
   * We will fill these before allocating new space.
   */
  update_unused_elements ();
  
  /**
   * Now refine the flagged elements.  This will
   * take up some space, maybe more than what was freed.
   */
  refine_elements();

  /**
   * Find any nodes that were orphaned by the coarsening
   * step and not reconnected during the refinement step.
   */
  update_unused_nodes();
  
  /**
   * We need contiguous elements, so
   * better trim the vector in case
   * we didn't fill up the coarsened
   * space.
   */
  mesh.trim_unused_elements (unused_elements);
  
  assert (unused_elements.empty());
  
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {
      assert (mesh.elem(e) != NULL);
      //mesh.elem(e)->set_refinement_flag(Elem::DO_NOTHING);
    }
  
  /**
   * Finally, the new mesh needs to be prepared for use
   */
  mesh.prepare_for_use ();
}



bool MeshRefinement::make_coarsening_compatible(const bool maintain_level_one)
{
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
   * sense to coarsen them we must de-set their coarsen flags if
   * they are set.
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {      
      assert (mesh.elem(e) != NULL);
      max_level = std::max(max_level, mesh.elem(e)->level());
      
      if (mesh.elem(e)->active()      &&
	  (mesh.elem(e)->level() == 0) &&
	  (mesh.elem(e)->refinement_flag() == Elem::COARSEN))
	mesh.elem(e)->set_refinement_flag(Elem::DO_NOTHING);
    }

  /**
   * if there are no refined elements then
   * there is no work for us to do
   */
  if (max_level == 0) return compatible_with_refinement; 

  
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
	  
	  active_elem_iterator       el    (mesh.elements_begin());
	  const active_elem_iterator end_el(mesh.elements_end());
	  
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
		      else  // I have a neighbor and it is not active. That means it has children.
			{   // While it _may_ be possible to coarsen us if all the children of
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
      
      /**
       * Next we look at all of the inactive cells.
       * If there is an inactive cell with all of its children
       * wanting to be unrefined then the element is a candidate
       * for unrefinement.  If all the children don't
       * all want to be unrefined then ALL of them need to have their
       * unrefinement flags cleared. 
       */
      for (int level=(max_level); level >= 0; level--)
	for (unsigned int e=0; e<mesh.n_elem(); e++)
	  if (!mesh.elem(e)->active() &&
	      (mesh.elem(e)->level() == static_cast<unsigned int>(level)))
	    {
	      /**
	       * right now the element hasn't been disqualified
	       * as a candidate for unrefinement
	       */
	      bool is_a_candidate = true;
	      
	      for (unsigned int c=0; c<mesh.elem(e)->n_children(); c++)
		if ((mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN) ||
		    !mesh.elem(e)->child(c)->active() )
		  is_a_candidate = false;
	      
	      if (!is_a_candidate)
		{
		  mesh.elem(e)->set_refinement_flag(Elem::DO_NOTHING);
		  
		  for (unsigned int c=0; c<mesh.elem(e)->n_children(); c++)
		    if (mesh.elem(e)->child(c)->refinement_flag() == Elem::COARSEN)
		      {
			level_one_satisfied = false;
			mesh.elem(e)->child(c)->set_refinement_flag(Elem::DO_NOTHING);
		      }
		}
	    }
      
      if (!level_one_satisfied) goto repeat;
      
    } // end if (maintian_level_one)

  
  /**
   * If all the children of a parent are set to be coarsened
   * then flag the parent so that they can kill thier kids...
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (!mesh.elem(e)->active())
      {
	/**
	 * Presume all the children are flagged for coarsening
	 * and then look for a contradiction
	 */
	bool all_children_flagged_for_coarsening = true;
	
	for (unsigned int c=0; c<mesh.elem(e)->n_children(); c++)
	  if (mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN)
	    all_children_flagged_for_coarsening = false;
	
	if (all_children_flagged_for_coarsening)
	  mesh.elem(e)->set_refinement_flag(Elem::COARSEN);
      }
	
  return compatible_with_refinement;
}



bool MeshRefinement::make_refinement_compatible(const bool maintain_level_one)
{
  // Unless we encounter a specific situation level-one
  // will be satisfied after executing this loop just once
  bool level_one_satisfied = true;

  // Unless we encounter a specific situation we will be compatible
  // with any selected coarsening flags
  bool compatible_with_coarsening = true;

#ifdef DEBUG
  
  // It only makes sense to refine active elements
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->refinement_flag() == Elem::REFINE)
      assert (mesh.elem(e)->active());
  
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
	  
	  for (unsigned int e=0; e<mesh.n_elem(); e++)
	    if (mesh.elem(e)->active() &&                
		(mesh.elem(e)->refinement_flag() == Elem::REFINE))  // If the element is active and the
	                                                            // refinement flag is set    
	      {
		const unsigned int my_level = mesh.elem(e)->level();
	    
		for (unsigned int side=0; side<mesh.elem(e)->n_sides(); side++)
		  if (mesh.elem(e)->neighbor(side) != NULL)     // I have a neighbor
		    if (mesh.elem(e)->neighbor(side)->active()) // and it is active
		      {
			Elem* neighbor = mesh.elem(e)->neighbor(side);
			
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

  
  return compatible_with_coarsening;
}



void MeshRefinement::coarsen_elements ()
{

  unsigned int max_level=0;

  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {
      assert (mesh.elem(e) != NULL);
      max_level = std::max(max_level, mesh.elem(e)->level());
    }
  
  if (max_level == 0) return;

  
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->refinement_flag() == Elem::COARSEN)
      {
	if (mesh.elem(e)->active())
	  {
	    assert (mesh.elem(e)->level() != 0);

	    // Remove this element from any neighbor
	    // lists that point to it.
	    mesh.elem(e)->nullify_neighbors();
	    
	    delete mesh.elem(e);
	    
	    mesh.add_elem(NULL, e);
	  }
	else
	  {
	    mesh.elem(e)->coarsen();
	    assert (mesh.elem(e)->active());
	  }
      }
}



void MeshRefinement::refine_elements ()
{
  const unsigned int orig_n_elem = mesh.n_elem();
  
  for (unsigned int e=0; e<orig_n_elem; e++)
    if (mesh.elem(e) != NULL) // could be if deleted
      if (mesh.elem(e)->refinement_flag() == Elem::REFINE)
	mesh.elem(e)->refine(mesh);
}



void MeshRefinement::uniformly_refine (unsigned int n)
{
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      for (unsigned int e=0; e<mesh.n_elem(); e++)
	if (mesh.elem(e)->active())
	  mesh.elem(e)->set_refinement_flag(Elem::REFINE);

      /**
       * Refine all the elements we just flagged.
       */
      refine_elements();
      
  
      /**
       * Finally, the new mesh needs to be prepared for use
       */
      mesh.prepare_for_use ();
    }
}


#endif
