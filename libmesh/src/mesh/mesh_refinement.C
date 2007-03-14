// $Id: mesh_refinement.C,v 1.57 2007-03-14 22:02:11 roystgnr Exp $

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



// Local includes
#include "libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef ENABLE_AMR

#include "error_vector.h"
#include "mesh_base.h"
#include "mesh_refinement.h"
#include "elem.h"
#include "libmesh_logging.h"
#include "boundary_info.h"



//-----------------------------------------------------------------
// Mesh refinement methods
MeshRefinement::MeshRefinement (MeshBase& m) :
  _mesh(m),
  _use_member_parameters(false),
  _coarsen_by_parents(false),
  _refine_fraction(0.3),
  _coarsen_fraction(0.0),
  _max_h_level(libMesh::invalid_uint),
  _coarsen_threshold(10),
  _nelem_target(0),
  _absolute_global_tolerance(0.0)
{
}



MeshRefinement::~MeshRefinement ()
{
  this->clear();  
}



void MeshRefinement::clear ()
{
  _new_nodes_map.clear();
}



Node* MeshRefinement::add_point (const Point& p,
                                 const unsigned int key,
                                 const Real tol)
{
  START_LOG("add_point()", "MeshRefinement");

  // Look for the key in the multimap  
  std::pair<map_type::iterator, map_type::iterator>
    pos = _new_nodes_map.equal_range(key);
  
      
  while (pos.first != pos.second) 
    if (p.absolute_fuzzy_equals(*(pos.first->second), tol))
    {
      STOP_LOG("add_point()", "MeshRefinement");
      return pos.first->second;      
    }
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
  STOP_LOG("add_point()", "MeshRefinement");
  return node;
}



Elem* MeshRefinement::add_elem (Elem* elem)
{
  assert (elem != NULL);

  
//   // If the unused_elements has any iterators from
//   // old elements, take the first one
//   if (!_unused_elements.empty())
//     {
//       std::vector<Elem*>::iterator it = _unused_elements.front();

//       *it = elem;

//       _unused_elements.pop_front();
//     }

//   // Otherwise, use the conventional add method
//   else
//     {
//       _mesh.add_elem (elem);
//     }

  // The _unused_elements optimization has been turned off.
  _mesh.add_elem (elem);
  
  return elem;
}



void MeshRefinement::create_parent_error_vector
  (const ErrorVector& error_per_cell,
   ErrorVector& error_per_parent,
   Real& parent_error_min,
   Real& parent_error_max)
{
  error_per_parent.clear();
  error_per_parent.resize(error_per_cell.size(), -1.);

  parent_error_min = std::numeric_limits<double>::max();
  parent_error_max = 0.;

  // We need to loop over all active elements to find minimum error
  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem   = *elem_it;
      Elem* parent = elem->parent();

      // Calculate summed errors on parent cells
      if (parent)
        {
          const unsigned int parentid  = parent->id();
          assert (parentid < error_per_parent.size());

          // If we haven't already calculated the parent's total
          // error, do so now
          if (error_per_parent[parentid] == -1.)
            {
              // The error estimator might have already given us an
              // estimate on the coarsenable parent elements
              if (error_per_cell[parentid])
                {
                  error_per_parent[parentid] = error_per_cell[parentid];
                  continue;
                }
 
              ErrorVectorReal parent_error = 0.;
              for (unsigned int n = 0; n != parent->n_children(); ++n)
                {
                  Elem* child = parent->child(n);
                  const unsigned int childid = child->id();

                  // If the parent has grandchildren we won't be able
                  // to coarsen it, so forget it
                  if (!child->active())
                    {
                      parent_error = -1.;
                      break;
                    }

                  // We take the square root of the sum of the
                  // squares, so errors that are Hilbert norms
                  // remain Hilbert norms
                  parent_error += (error_per_cell[childid] *
                                   error_per_cell[childid]);
                }

              // If this element is uncoarsenable, just skip it
              if (parent_error < 0)
                continue;

              parent_error = sqrt(parent_error);
              error_per_parent[parentid] = parent_error;

              parent_error_min = std::min (parent_error_min,
                                           parent_error);
              parent_error_max = std::max (parent_error_max,
                                           parent_error);
            }
        }
    }
}



void MeshRefinement::update_nodes_map ()
{
  START_LOG("update_nodes_map()", "MeshRefinement");

  _new_nodes_map.clear();

  MeshBase::node_iterator       it  = _mesh.nodes_begin();
  const MeshBase::node_iterator end = _mesh.nodes_end();

  for (; it != end; ++it)
    {
      Node* node = *it;

      // Add the node to the map.
      _new_nodes_map.insert(std::make_pair(node->key(), node));
    }

  STOP_LOG("update_nodes_map()", "MeshRefinement");
} 



bool MeshRefinement::test_level_one (bool assert_pass)
{
  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem *elem = *elem_it;

      for (unsigned int n=0; n<elem->n_neighbors(); n++)
        {
          Elem *neighbor = elem->neighbor(n);

          if (!neighbor || !neighbor->active())
            continue;

          if ((neighbor->level() + 1 < elem->level()) ||
              (neighbor->p_level() + 1 < elem->p_level()) ||
              (neighbor->p_level() > elem->p_level() + 1))
            {
              // We didn't pass the level one test, so assert that 
              // we're allowed not to
              assert(!assert_pass);
              return false;
            }
        }
    }
  return true;
}



bool MeshRefinement::test_unflagged (bool assert_pass)
{
  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem *elem = *elem_it;

      if (elem->refinement_flag() == Elem::REFINE ||
          elem->refinement_flag() == Elem::COARSEN ||
          elem->p_refinement_flag() == Elem::REFINE ||
          elem->p_refinement_flag() == Elem::COARSEN)
        {
          // We didn't pass the "elements are unflagged" test,
          // so assert that we're allowed not to
          assert(!assert_pass);
          return false;
        }
    }
  return true;
}



bool MeshRefinement::refine_and_coarsen_elements (const bool maintain_level_one)
{
  //assert (_mesh.mesh_dimension() != 1);

  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (maintain_level_one)
    assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem *elem = *elem_it;

      // Set refinement flag to INACTIVE if the
      // element isn't active
      if ( !elem->active())
        {
	  elem->set_refinement_flag(Elem::INACTIVE);
	  elem->set_p_refinement_flag(Elem::INACTIVE);
        }

      // This might be left over from the last step
      if (elem->refinement_flag() == Elem::JUST_REFINED)
	elem->set_refinement_flag(Elem::DO_NOTHING);
    }

  
  // Repeat until coarsening & refinement flags jive
  bool satisfied = false;
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

  // First coarsen the flagged elements.
  const bool coarsening_changed_mesh =
    this->_coarsen_elements ();

  if (maintain_level_one)
    assert(test_level_one(true));
  assert(this->make_coarsening_compatible(maintain_level_one));
  assert(this->make_refinement_compatible(maintain_level_one));
// FIXME: This won't pass unless we add a redundant find_neighbors()
// call or replace find_neighbors() with on-the-fly neighbor updating
// assert(!this->eliminate_unrefined_patches());

  // We can't contract the mesh ourselves anymore - a System might
  // need to restrict old coefficient vectors first
  // _mesh.contract();

  // Now refine the flagged elements.  This will
  // take up some space, maybe more than what was freed.
  const bool refining_changed_mesh =
    this->_refine_elements();
  
  if (maintain_level_one)
    assert(test_level_one(true));
  assert(test_unflagged(true));
  assert(this->make_coarsening_compatible(maintain_level_one));
  assert(this->make_refinement_compatible(maintain_level_one));
// FIXME: This won't pass unless we add a redundant find_neighbors()
// call or replace find_neighbors() with on-the-fly neighbor updating
// assert(!this->eliminate_unrefined_patches());

  // Finally, the new mesh needs to be prepared for use
  if (coarsening_changed_mesh || refining_changed_mesh)
    {
      _mesh.prepare_for_use ();
      
      return true;
    }

  // Otherwise there was no change in the mesh,
  // let the user know.  Also, there is no need
  // to prepare the mesh for use since it did not change.
  return false;
    
}







bool MeshRefinement::coarsen_elements (const bool maintain_level_one)
{
  //assert (_mesh.mesh_dimension() != 1);
  
  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (maintain_level_one)
    assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem* elem = *elem_it;

      // Set refinement flag to INACTIVE if the
      // element isn't active
      if ( !elem->active())
        {
	  elem->set_refinement_flag(Elem::INACTIVE);
	  elem->set_p_refinement_flag(Elem::INACTIVE);
        }

      // This might be left over from the last step
      if (elem->refinement_flag() == Elem::JUST_REFINED)
	elem->set_refinement_flag(Elem::DO_NOTHING);
    }

  // Repeat until the flags form a conforming mesh.
  bool satisfied = false;
  while (!satisfied)
    {
      const bool coarsening_satisfied =
	this->make_coarsening_compatible(maintain_level_one);

      const bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();// &&
// 	!this->limit_level_mismatch_at_node(1);
      
      satisfied = (coarsening_satisfied &&
		   smoothing_satisfied);
    }

  // Coarsen the flagged elements.
  const bool mesh_changed = 
    this->_coarsen_elements ();

  if (maintain_level_one)
    assert(test_level_one(true));
  assert(this->make_coarsening_compatible(maintain_level_one));
// FIXME: This won't pass unless we add a redundant find_neighbors()
// call or replace find_neighbors() with on-the-fly neighbor updating
// assert(!this->eliminate_unrefined_patches());
    
  // We can't contract the mesh ourselves anymore - a System might
  // need to restrict old coefficient vectors first
  // _mesh.contract();

  // Finally, the new mesh may need to be prepared for use
  if (mesh_changed)
    _mesh.prepare_for_use ();

  return mesh_changed;
}







bool MeshRefinement::refine_elements (const bool maintain_level_one)
{
  //assert (_mesh.mesh_dimension() != 1);

  if (maintain_level_one)
    assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem *elem = *elem_it;

      // Set refinement flag to INACTIVE if the
      // element isn't active
      if ( !elem->active())
        {
	  elem->set_refinement_flag(Elem::INACTIVE);
	  elem->set_p_refinement_flag(Elem::INACTIVE);
        }

      // This might be left over from the last step
      if (elem->refinement_flag() == Elem::JUST_REFINED)
	elem->set_refinement_flag(Elem::DO_NOTHING);
    }

  
  // Repeat until coarsening & refinement flags jive
  bool satisfied = false;
  while (!satisfied)
    {
      const bool refinement_satisfied =
	this->make_refinement_compatible(maintain_level_one);

      const bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();// &&
// 	!this->limit_level_mismatch_at_node(1);
      
      satisfied = (refinement_satisfied &&
		   smoothing_satisfied);
    }
  
  // Now refine the flagged elements.  This will
  // take up some space, maybe more than what was freed.
  const bool mesh_changed = 
    this->_refine_elements();

  if (maintain_level_one)
    assert(test_level_one(true));
  assert(this->make_refinement_compatible(maintain_level_one));
// FIXME: This won't pass unless we add a redundant find_neighbors()
// call or replace find_neighbors() with on-the-fly neighbor updating
// assert(!this->eliminate_unrefined_patches());
    
  // Finally, the new mesh needs to be prepared for use
  if (mesh_changed)
    _mesh.prepare_for_use ();

  return mesh_changed;
}






bool MeshRefinement::make_coarsening_compatible(const bool maintain_level_one)
{
  START_LOG ("make_coarsening_compatible()", "MeshRefinement");
  
  
  // Unless we encounter a specific situation level-one
  // will be satisfied after executing this loop just once   
  bool level_one_satisfied = true;

  
  // Unless we encounter a specific situation we will be compatible
  // with any selected refinement flags
  bool compatible_with_refinement = true;

  
  // find the maximum h and p levels in the mesh   
  unsigned int max_level = 0;
  unsigned int max_p_level = 0;
    
  
  // First we look at all the active level-0 elements.  Since it doesn't make
  // sense to coarsen them we must un-set their coarsen flags if
  // they are set.   
  MeshBase::element_iterator       el     = _mesh.active_elements_begin();
  const MeshBase::element_iterator end_el = _mesh.active_elements_end(); 

  for (; el != end_el; ++el)
    {      
      Elem *elem = *el;
      max_level = std::max(max_level, elem->level());
      max_p_level =
        std::max(max_p_level,
		 static_cast<unsigned int>(elem->p_level()));
      
      if ((elem->level() == 0) &&
	  (elem->refinement_flag() == Elem::COARSEN))
	elem->set_refinement_flag(Elem::DO_NOTHING);

      if ((elem->p_level() == 0) &&
	  (elem->p_refinement_flag() == Elem::COARSEN))
	elem->set_p_refinement_flag(Elem::DO_NOTHING);
    }
  
  // if there are no refined elements then
  // there is no work for us to do   
  if (max_level == 0 && max_p_level == 0)
    {
      STOP_LOG ("make_coarsening_compatible()", "MeshRefinement");
      return compatible_with_refinement;
    }

  
  
  // Loop over all the active elements.  If an element is marked
  // for coarsening we better check its neighbors.  If ANY of these neighbors
  // are marked for refinement AND are at the same level then there is a
  // conflict.  By convention refinement wins, so we un-mark the element for
  // coarsening.  Level-one would be violated in this case so we need to re-run
  // the loop.   
  if (maintain_level_one)
    {
      
    repeat:
      level_one_satisfied = true;
  
      do
	{
	  level_one_satisfied = true;

	  MeshBase::element_iterator       el     = _mesh.active_elements_begin();
	  const MeshBase::element_iterator end_el = _mesh.active_elements_end(); 

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
                              break;
			    }
			}
		      else // I have a neighbor and it is not active. That means it has children.
			{  // While it _may_ be possible to coarsen us if all the children of
			   // that element want to be coarsened, it is impossible to know at this
			   // stage.  Forget about it for the moment...  This can be handled in
			   // two steps.
			  elem->set_refinement_flag(Elem::DO_NOTHING);
			  level_one_satisfied = false;
                          break;
			}
		}
	      if (elem->p_refinement_flag() == Elem::COARSEN) // If
                // the element is active and the order reduction flag is set
		{
		  const unsigned int my_p_level = elem->p_level();

		  for (unsigned int n=0; n<elem->n_neighbors(); n++)
		    if (elem->neighbor(n) != NULL)     // I have a neighbor
		      if (elem->neighbor(n)->active()) // and it is active
			{
			  const Elem* neighbor = elem->neighbor(n);

                          if ((neighbor->p_level() > my_p_level &&
                               neighbor->p_refinement_flag() != Elem::COARSEN)
                              || (neighbor->p_level() == my_p_level &&
                               neighbor->p_refinement_flag() == Elem::REFINE))
			    {
			      elem->set_p_refinement_flag(Elem::DO_NOTHING);
			      level_one_satisfied = false;
                              break;
			    }
			}
		      else // I have a neighbor and it is not active.
			{  // We need to find which of its children
                           // have me as a neighbor, and maintain
                           // level one p compatibility with them.
                           // Because we currently have level one h
                           // compatibility, we don't need to check
                           // grandchildren
			  const Elem* neighbor = elem->neighbor(n);

                           assert(neighbor->has_children());
	                   for (unsigned int c=0; c!=neighbor->n_children(); c++)
                             {
                               Elem *subneighbor = neighbor->child(c);
                               if (subneighbor->active() &&
                                   subneighbor->is_neighbor(elem))
                                 if ((subneighbor->p_level() > my_p_level &&
                                     subneighbor->p_refinement_flag() != Elem::COARSEN)
                                     || (subneighbor->p_level() == my_p_level &&
                                     subneighbor->p_refinement_flag() == Elem::REFINE))
			           {
			             elem->set_p_refinement_flag(Elem::DO_NOTHING);
			             level_one_satisfied = false;
                                     break;
			           }
                             }
                           if (!level_one_satisfied)
                             break;
			}
		}
		  
	    }      
	}
      while (!level_one_satisfied);
      
    } // end if (maintain_level_one)
  
  
  // Next we look at all of the ancestor cells.
  // If there is a parent cell with all of its children
  // wanting to be unrefined then the element is a candidate
  // for unrefinement.  If all the children don't
  // all want to be unrefined then ALL of them need to have their
  // unrefinement flags cleared.    
  for (int level=(max_level); level >= 0; level--)
    for (unsigned int e=0; e<_mesh.n_elem(); e++)
    {
      if (_mesh.elem(e)->ancestor() &&
	  (_mesh.elem(e)->level() == static_cast<unsigned int>(level)))
	{
	  
	  // right now the element hasn't been disqualified
	  // as a candidate for unrefinement	   
	  bool is_a_candidate = true;
	      
	  for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
	    if ((_mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN) ||
		!_mesh.elem(e)->child(c)->active() )
	      is_a_candidate = false;
	      
	  if (!is_a_candidate)
	    {
	      _mesh.elem(e)->set_refinement_flag(Elem::INACTIVE);
		  
	      for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
		if (_mesh.elem(e)->child(c)->refinement_flag() == Elem::COARSEN)
		  {
		    level_one_satisfied = false;
		    _mesh.elem(e)->child(c)->set_refinement_flag(Elem::DO_NOTHING);
		  }
	    }
	}
     }
     
  if (!level_one_satisfied && maintain_level_one) goto repeat;
  
  
  // If all the children of a parent are set to be coarsened
  // then flag the parent so that they can kill thier kids...   
  for (unsigned int e=0; e<_mesh.n_elem(); e++)
    if (_mesh.elem(e)->ancestor())
      {
	
	// Presume all the children are flagged for coarsening
	// and then look for a contradiction	 
	bool all_children_flagged_for_coarsening = true;
	
	for (unsigned int c=0; c<_mesh.elem(e)->n_children(); c++)
	  if (_mesh.elem(e)->child(c)->refinement_flag() != Elem::COARSEN)
	    all_children_flagged_for_coarsening = false;
	
	if (all_children_flagged_for_coarsening)
	  _mesh.elem(e)->set_refinement_flag(Elem::COARSEN_INACTIVE);
        else
	  _mesh.elem(e)->set_refinement_flag(Elem::INACTIVE);
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

  // This loop enforces the level-1 rule.  We should only
  // execute it if the user indeed wants level-1 satisfied!   
  if (maintain_level_one)
    {  
      do
	{
	  level_one_satisfied = true;
	  
	  MeshBase::element_iterator       el     = _mesh.active_elements_begin();
	  const MeshBase::element_iterator end_el = _mesh.active_elements_end(); 

	  for (; el != end_el; ++el)
            {
            Elem *elem = *el;
	    if (elem->refinement_flag() == Elem::REFINE)  // If the element is active and the
                                                          // h refinement flag is set    
	      {
		const unsigned int my_level = elem->level();
	    
		for (unsigned int side=0; side != elem->n_sides(); side++)
                  {
		    Elem* neighbor = elem->neighbor(side);
		    if (neighbor != NULL &&     // I have a neighbor
		        neighbor->active()) // and it is active
		      {
			
			
			// Case 1:  The neighbor is at the same level I am.
			//        1a: The neighbor will be refined       -> NO PROBLEM
			//        1b: The neighbor won't be refined      -> NO PROBLEM
			//        1c: The neighbor wants to be coarsened -> PROBLEM			 
			if (neighbor->level() == my_level)
			  {
			    if (neighbor->refinement_flag() == Elem::COARSEN)
			      {
				neighbor->set_refinement_flag(Elem::DO_NOTHING);
                                if (neighbor->parent())
                                  neighbor->parent()->set_refinement_flag(Elem::INACTIVE);
				compatible_with_coarsening = false;
				level_one_satisfied = false;
			      }
			  }
			
			
			// Case 2: The neighbor is one level lower than I am.
			//         The neighbor thus MUST be refined to satisfy
			//         the level-one rule, regardless of whether it
			//         was originally flagged for refinement. If it
			//         wasn't flagged already we need to repeat
			//         this process.			 
			else if ((neighbor->level()+1) == my_level)
			  {
			    if (neighbor->refinement_flag() != Elem::REFINE)
			      {
				neighbor->set_refinement_flag(Elem::REFINE);
                                if (neighbor->parent())
				  neighbor->parent()->set_refinement_flag(Elem::INACTIVE);
				compatible_with_coarsening = false;
				level_one_satisfied = false; 
			      }
			  }
#ifdef DEBUG		
			
			// Sanity check. We should never get into a
			// case when our neighbot is more than one
			// level away.			 
			else if ((neighbor->level()+1) < my_level)
			  {
			    error();
			  }
			
			
			// Note that the only other possibility is that the
			// neighbor is already refined, in which case it isn't
			// active and we should never get here.			 
			else
			  {
			    error();
			  }
#endif		
		      } 
                  }
	      }
	    if (elem->p_refinement_flag() == Elem::REFINE)  // If the element is active and the
                                                            // p refinement flag is set    
	      {
		const unsigned int my_p_level = elem->p_level();

		for (unsigned int side=0; side != elem->n_sides(); side++)
                  {
                    Elem *neighbor = elem->neighbor(side);
		    if (neighbor != NULL)     // I have a neighbor
		      if (neighbor->active()) // and it is active
		        {
                          if (neighbor->p_level() < my_p_level &&
                               neighbor->p_refinement_flag() != Elem::REFINE)
			    {
			      neighbor->set_p_refinement_flag(Elem::REFINE);
			      level_one_satisfied = false;
			      compatible_with_coarsening = false;
			    }
                           if (neighbor->p_level() == my_p_level &&
                              neighbor->p_refinement_flag() == Elem::COARSEN)
			    {
			      neighbor->set_p_refinement_flag(Elem::DO_NOTHING);
			      level_one_satisfied = false;
			      compatible_with_coarsening = false;
			    }
		        } 
		      else // I have an inactive neighbor
		        {
                           assert(neighbor->has_children());
	                   for (unsigned int c=0; c!=neighbor->n_children(); c++)
                             {
                               Elem *subneighbor = neighbor->child(c);
                               if (subneighbor->active() &&
                                   subneighbor->is_neighbor(elem))
                                 if (subneighbor->p_level() < my_p_level &&
                                     subneighbor->p_refinement_flag() != Elem::REFINE)
			           {
                                     // We should already be level one
                                     // compatible
                                     assert(subneighbor->p_level() + 2u > 
                                            my_p_level);
			             subneighbor->set_p_refinement_flag(Elem::REFINE);
			             level_one_satisfied = false;
			             compatible_with_coarsening = false;
			           }
                                 if (subneighbor->p_level() == my_p_level &&
                                     subneighbor->p_refinement_flag() == Elem::COARSEN)
			           {
			             subneighbor->set_p_refinement_flag(Elem::DO_NOTHING);
			             level_one_satisfied = false;
			             compatible_with_coarsening = false;
			           }
                             }
		        } 
                  }
	      }
            }
	}
      
      while (!level_one_satisfied);
    } // end if (maintain_level_one)

  
  STOP_LOG ("make_refinement_compatible()", "MeshRefinement");
  
  return compatible_with_coarsening;
}




bool MeshRefinement::_coarsen_elements ()
{
  START_LOG ("_coarsen_elements()", "MeshRefinement");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;
  
  // Clear the unused_elements data structure.
  // The elements have been packed since it was built,
  // so there are _no_ unused elements.  We cannot trust
  // any iterators currently in this data structure.
  // _unused_elements.clear();

  MeshBase::element_iterator       it  = _mesh.elements_begin();
  const MeshBase::element_iterator end = _mesh.elements_end();

  // Loop over the elements.   
  for ( ; it != end; ++it)
    {
      Elem* elem = *it;

      // Not necessary when using elem_iterator
      // assert (elem != NULL);
      
      // active elements flagged for coarsening will
      // no longer be deleted until MeshRefinement::contract()
      if (elem->refinement_flag() == Elem::COARSEN)
	{
	  // Huh?  no level-0 element should be active
	  // and flagged for coarsening.
	  assert (elem->level() != 0);
	  
	  // Remove this element from any neighbor
	  // lists that point to it.
	  elem->nullify_neighbors();

	  // Remove any boundary information associated
	  // with this element
	  _mesh.boundary_info->remove (elem);
	      
	  // Add this iterator to the _unused_elements
	  // data structure so we might fill it.
	  // The _unused_elements optimization is currently off.
	  // _unused_elements.push_back (it);

          // About to call parent()->active!

	  // Is this element's parent inactive??
	  if (!elem->parent()->active())
	    {
	      std::cerr << "Element being deleted has inactive parent." << std::endl;
	    }

	  // Don't delete the element until
	  // MeshRefinement::contract()
	  // _mesh.delete_elem(elem);

	  // the mesh has certainly changed
	  mesh_changed = true;
        }
      
      // inactive elements flagged for coarsening
      // will become active
      else if (elem->refinement_flag() == Elem::COARSEN_INACTIVE)
	{
	  elem->coarsen();
	  assert (elem->active());

	  // the mesh has certainly changed
	  mesh_changed = true;
	}
      if (elem->p_refinement_flag() == Elem::COARSEN)
        {
          if (elem->p_level() > 0)
            {
              elem->set_p_refinement_flag(Elem::JUST_COARSENED);
              elem->set_p_level(elem->p_level() - 1);
	      mesh_changed = true;
            }
          else
            {
              elem->set_p_refinement_flag(Elem::DO_NOTHING);
            }
        }
    }  
  
  STOP_LOG ("_coarsen_elements()", "MeshRefinement");

  return mesh_changed;
}







bool MeshRefinement::_refine_elements ()
{
  // Update the _new_nodes_map so that elements can
  // find nodes to connect to.
  this->update_nodes_map ();

  START_LOG ("_refine_elements()", "MeshRefinement");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  // Get the original number of elements.
  const unsigned int orig_n_elem = _mesh.n_elem();

  // Construct a local vector of Elem* which have been
  // previously marked for refinement.  We reserve enough
  // space to allow for every element to be refined.
  std::vector<Elem*> local_copy_of_elements;
  local_copy_of_elements.reserve(orig_n_elem);

  // Iterate over the elements, looking for elements
  // flagged for refinement.
  MeshBase::element_iterator       it  = _mesh.elements_begin();
  const MeshBase::element_iterator end = _mesh.elements_end();

  for (; it != end; ++it)
    {
      Elem* elem = *it;
      if (elem->refinement_flag() == Elem::REFINE)
	local_copy_of_elements.push_back(elem);
      if (elem->p_refinement_flag() == Elem::REFINE &&
          elem->active())
        {
	  elem->set_p_level(elem->p_level()+1);
	  elem->set_p_refinement_flag(Elem::JUST_REFINED);
          mesh_changed = true;
        }
    }

  // The mesh will change if there are elements to refine
  if(!(local_copy_of_elements.empty()))
    mesh_changed = true;
  
  // Now iterate over the local copy and refine each one.
  // This may resize the mesh's internal container and invalidate
  // any existing iterators.
  for (unsigned int e=0; e<local_copy_of_elements.size(); ++e)
    local_copy_of_elements[e]->refine(*this);
  
  // Clear the _new_nodes_map and _unused_elements data structures.
  this->clear();
  
  STOP_LOG ("_refine_elements()", "MeshRefinement");

  return mesh_changed;
}



void MeshRefinement::uniformly_p_refine (unsigned int n)
{
  // Refine n times
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // P refine all the active elements
      MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
        {
	  (*elem_it)->set_p_level((*elem_it)->p_level()+1);
	  (*elem_it)->set_p_refinement_flag(Elem::JUST_REFINED);
        }
    }
}



void MeshRefinement::uniformly_p_coarsen (unsigned int n)
{
  // Coarsen p times
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // P coarsen all the active elements
      MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
        {
          if ((*elem_it)->p_level() > 0)
            {
	      (*elem_it)->set_p_level((*elem_it)->p_level()-1);
	      (*elem_it)->set_p_refinement_flag(Elem::JUST_COARSENED);
            }
        }
    }
}



void MeshRefinement::uniformly_refine (unsigned int n)
{
  // Refine n times
  // FIXME - this won't work if n>1 and the mesh 
  // has already been attached to an equation system
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // Clean up the refinement flags
      this->clean_refinement_flags();
      
      // Flag all the active elements for refinement.       
      MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_refinement_flag(Elem::REFINE);

      // Refine all the elements we just flagged.
      this->_refine_elements();
    }
  
  // Finally, the new mesh needs to be prepared for use
  _mesh.prepare_for_use ();
}



void MeshRefinement::uniformly_coarsen (unsigned int n)
{
  // Coarsen n times
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // Clean up the refinement flags
      this->clean_refinement_flags();
      
      // Flag all the active elements for coarsening
      MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
        {
	  (*elem_it)->set_refinement_flag(Elem::COARSEN);
          if ((*elem_it)->parent()) 
	    (*elem_it)->parent()->set_refinement_flag(Elem::COARSEN_INACTIVE);
        }

      // Coarsen all the elements we just flagged.
      this->_coarsen_elements();
    }
    
  
  // Finally, the new mesh needs to be prepared for use
  _mesh.prepare_for_use ();
}


#endif
