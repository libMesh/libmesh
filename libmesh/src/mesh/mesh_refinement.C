// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <limits>

// Local includes
#include "libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef ENABLE_AMR

#include "boundary_info.h"
#include "elem.h"
#include "error_vector.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "mesh_refinement.h"
#include "parallel.h"



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
  _absolute_global_tolerance(0.0),
  _face_level_mismatch_limit(1),
  _edge_level_mismatch_limit(0),
  _node_level_mismatch_limit(0)
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
                                 const unsigned int processor_id,
                                 const Real tol)
{
  START_LOG("add_point()", "MeshRefinement");

  const unsigned int key = this->point_key(p);

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

  // Add the node, with a default id and the requested
  // processor_id
  Node* node = _mesh.add_point (p, DofObject::invalid_id,
                                processor_id);

  assert (node != NULL);

  // Add the node to the map.  In the case of the
  // std::multimap use pos.first as a hint for where to put it
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
  _new_nodes_map.insert(std::make_pair(key, node));
#else
  _new_nodes_map.insert(pos.first, std::make_pair(key, node));
#endif			    

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
  // This function must be run on all processors at once
  parallel_only();

  error_per_parent.clear();
  error_per_parent.resize(error_per_cell.size(), 0);

  // The parent's error is defined as the square root of the
  // sum of the children's errors squared, so errors that are
  // Hilbert norms remain Hilbert norms.
  //
  // Because the children may be on different processors, we
  // calculate local contributions to the parents' errors squared
  // first, then sum across processors and take the square roots
  // second.
  MeshBase::element_iterator       elem_it  = _mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem   = *elem_it;
      Elem* parent = elem->parent();

      // Calculate each contribution to parent cells
      if (parent)
        {
          const unsigned int parentid  = parent->id();
          assert (parentid < error_per_parent.size());

          // If the parent has grandchildren we won't be able
          // to coarsen it, so forget it.
          // On a parallel mesh, all the parent's children will be
          // local or ghost elements - it's grandchildren may be
	  // remote elements, but that's fine; we only need to know if
	  // they exist.
	  bool parent_has_grandchildren = false;
          for (unsigned int n = 0; n != parent->n_children(); ++n)
            {
              Elem* child = parent->child(n);

              if (!child->active())
                {
                  parent_has_grandchildren = true;
                  break;
                }
            }

          if (!parent_has_grandchildren)
            error_per_parent[parentid] += (error_per_cell[elem->id()] *
                                           error_per_cell[elem->id()]);
        }
    }

  // Sum the vector across all processors
  Parallel::sum(static_cast<std::vector<ErrorVectorReal>&>(error_per_parent));

  // Calculate the min and max as we loop
  parent_error_min = std::numeric_limits<double>::max();
  parent_error_max = 0.;

  for (unsigned int i = 0; i != error_per_parent.size(); ++i)
    {
      // The error estimator might have already given us an
      // estimate on the coarsenable parent elements; if so then
      // we want to retain that estimate
      if (error_per_cell[i])
        {
          error_per_parent[i] = error_per_cell[i];
          continue;
        }
 
      // If this element isn't a coarsenable parent with error, we
      // have nothing to do.
      if (!error_per_parent[i])
        continue;

      error_per_parent[i] = std::sqrt(error_per_parent[i]);

      parent_error_min = std::min (parent_error_min,
                                   error_per_parent[i]);
      parent_error_max = std::max (parent_error_max,
                                   error_per_parent[i]);
    }
}



void MeshRefinement::update_nodes_map ()
{
  // This function must be run on all processors at once
  // for non-serial meshes
  if (!_mesh.is_serial())
    parallel_only();

  START_LOG("update_nodes_map()", "MeshRefinement");

  // Clear the old map
  _new_nodes_map.clear();

  // Cache a bounding box
  _lower_bound.clear();
  _lower_bound.resize(3, std::numeric_limits<Real>::max());
  _upper_bound.clear();
  _upper_bound.resize(3, -std::numeric_limits<Real>::max());

  MeshBase::node_iterator       it  = _mesh.nodes_begin();
  const MeshBase::node_iterator end = _mesh.nodes_end();

  for (; it != end; ++it)
    {
      Node* node = *it;

      // Expand the bounding box if necessary
      _lower_bound[0] = std::min(_lower_bound[0],
                                 (*node)(0));
      _lower_bound[1] = std::min(_lower_bound[1],
                                 (*node)(1));
      _lower_bound[2] = std::min(_lower_bound[2],
                                 (*node)(2));
      _upper_bound[0] = std::max(_upper_bound[0],
                                 (*node)(0));
      _upper_bound[1] = std::max(_upper_bound[1],
                                 (*node)(1));
      _upper_bound[2] = std::max(_upper_bound[2],
                                 (*node)(2));
    }

  // On a parallel mesh we might not yet have a full bounding box
  if (!_mesh.is_serial())
    {
      Parallel::min(_lower_bound);
      Parallel::max(_upper_bound);
    }

  // Populate the nodes map
  it  = _mesh.nodes_begin();

  for (; it != end; ++it)
    {
      Node* node = *it;

      // Add the node to the map.
      _new_nodes_map.insert(std::make_pair(this->point_key(*node), node));
    }

  STOP_LOG("update_nodes_map()", "MeshRefinement");
} 



bool MeshRefinement::test_level_one (bool assert_pass)
{
  // This function must be run on all processors at once
  parallel_only();

  MeshBase::element_iterator       elem_it  = _mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_local_elements_end();

  bool failure = false;

  for ( ; elem_it != elem_end && !failure; ++elem_it)
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
              failure = true;
              break;
            }
        }
    }

  // If any processor failed, we failed globally
  Parallel::max(failure);

  if (failure)
    {
      // We didn't pass the level one test, so assert that 
      // we're allowed not to
      assert(!assert_pass);
      return false;
    }
  return true;
}



bool MeshRefinement::test_unflagged (bool assert_pass)
{
  // This function must be run on all processors at once
  parallel_only();

  bool found_flag = false;

  // Search for local flags
  MeshBase::element_iterator       elem_it  = _mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_local_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Pointer to the element
      Elem *elem = *elem_it;

      if (elem->refinement_flag() == Elem::REFINE ||
          elem->refinement_flag() == Elem::COARSEN ||
          elem->p_refinement_flag() == Elem::REFINE ||
          elem->p_refinement_flag() == Elem::COARSEN)
        {
          found_flag = true;
          break;
        }
    }

  // If we found a flag on any processor, it counts
  Parallel::max(found_flag);

  if (found_flag)
    {
      // We didn't pass the "elements are unflagged" test,
      // so assert that we're allowed not to
      assert(!assert_pass);
      return false;
    }
  return true;
}



bool MeshRefinement::refine_and_coarsen_elements (const bool maintain_level_one)
{
  // This function must be run on all processors at once
  parallel_only();

  bool _maintain_level_one = maintain_level_one;

  // If the user used non-default parameters, let's warn that they're
  // deprecated
  if (!maintain_level_one)
    {
      deprecated();
    }
  else
    _maintain_level_one = _face_level_mismatch_limit;

  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (_maintain_level_one)
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

  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
  if (!_mesh.is_serial())
    this->make_flags_parallel_consistent();

  // Repeat until coarsening & refinement flags jive
  bool satisfied = false;
  do
    {
      const bool coarsening_satisfied =
	this->make_coarsening_compatible(maintain_level_one);

      const bool refinement_satisfied =
	this->make_refinement_compatible(maintain_level_one);

      bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();

      if (_edge_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied && 
          !this->limit_level_mismatch_at_edge (_edge_level_mismatch_limit);

      if (_node_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied &&
          !this->limit_level_mismatch_at_node (_node_level_mismatch_limit);

      // Parallel consistency has to come first, or coarsening
      // along processor boundaries might occasionally be falsely
      // prevented

      const bool parallel_consistent = _mesh.is_serial() ||
        this->make_flags_parallel_consistent();

      satisfied = (parallel_consistent &&
                   coarsening_satisfied &&
		   refinement_satisfied &&
		   smoothing_satisfied);
#ifdef DEBUG
      bool max_satisfied = satisfied,
           min_satisfied = satisfied;
      Parallel::max(max_satisfied);
      Parallel::min(min_satisfied);
      assert (satisfied == max_satisfied);
      assert (satisfied == min_satisfied);
#endif
    }
  while (!satisfied);

  // First coarsen the flagged elements.
  const bool coarsening_changed_mesh =
    this->_coarsen_elements ();

  if (_maintain_level_one)
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

  if (_maintain_level_one)
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
  // This function must be run on all processors at once
  parallel_only();

  bool _maintain_level_one = maintain_level_one;

  // If the user used non-default parameters, let's warn that they're
  // deprecated
  if (!maintain_level_one)
    {
      deprecated();
    }
  else
    _maintain_level_one = _face_level_mismatch_limit;

  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (_maintain_level_one)
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

  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
  if (!_mesh.is_serial())
    this->make_flags_parallel_consistent();

  // Repeat until the flags form a conforming mesh.
  bool satisfied = false;
  while (!satisfied)
    {
      const bool coarsening_satisfied =
	this->make_coarsening_compatible(maintain_level_one);

      bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();// &&

      if (_edge_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied &&
          !this->limit_level_mismatch_at_edge (_edge_level_mismatch_limit);

      if (_node_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied &&
          !this->limit_level_mismatch_at_node (_node_level_mismatch_limit);

      const bool parallel_consistent = _mesh.is_serial() ||
        this->make_flags_parallel_consistent();
      
      satisfied = (parallel_consistent &&
                   coarsening_satisfied &&
		   smoothing_satisfied);
#ifdef DEBUG
      bool max_satisfied = satisfied,
           min_satisfied = satisfied;
      Parallel::max(max_satisfied);
      Parallel::min(min_satisfied);
      assert (satisfied == max_satisfied);
      assert (satisfied == min_satisfied);
#endif
    }

  // Coarsen the flagged elements.
  const bool mesh_changed = 
    this->_coarsen_elements ();

  if (_maintain_level_one)
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
  // This function must be run on all processors at once
  parallel_only();

  bool _maintain_level_one = maintain_level_one;

  // If the user used non-default parameters, let's warn that they're
  // deprecated
  if (!maintain_level_one)
    {
      deprecated();
    }
  else
    _maintain_level_one = _face_level_mismatch_limit;

  if (_maintain_level_one)
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


  
  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
  if (!_mesh.is_serial())
    this->make_flags_parallel_consistent();

  // Repeat until coarsening & refinement flags jive
  bool satisfied = false;
  while (!satisfied)
    {
      const bool refinement_satisfied =
	this->make_refinement_compatible(maintain_level_one);

      bool smoothing_satisfied = 
 	!this->eliminate_unrefined_patches();// &&

      if (_edge_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied &&
          !this->limit_level_mismatch_at_edge (_edge_level_mismatch_limit);

      if (_node_level_mismatch_limit)
        smoothing_satisfied = smoothing_satisfied &&
          !this->limit_level_mismatch_at_node (_node_level_mismatch_limit);

      const bool parallel_consistent = _mesh.is_serial() ||
        this->make_flags_parallel_consistent();
      
      satisfied = (parallel_consistent &&
                   refinement_satisfied &&
		   smoothing_satisfied);
#ifdef DEBUG
      bool max_satisfied = satisfied,
           min_satisfied = satisfied;
      Parallel::max(max_satisfied);
      Parallel::min(min_satisfied);
      assert (satisfied == max_satisfied);
      assert (satisfied == min_satisfied);
#endif
    }
  
  // Now refine the flagged elements.  This will
  // take up some space, maybe more than what was freed.
  const bool mesh_changed = 
    this->_refine_elements();

  if (_maintain_level_one)
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






bool MeshRefinement::make_flags_parallel_consistent()
{
  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("make_flags_parallel_consistent()", "MeshRefinement");
  // We're consistent until we discover otherwise
  bool parallel_consistent = true;

  // Count the number of ghost elements we'll need to inquire about
  // from each other processor
  std::vector<unsigned int>
    ghost_elems_from_proc(libMesh::n_processors(), 0);

  const MeshBase::element_iterator end = _mesh.elements_end();

  for (MeshBase::element_iterator it = _mesh.elements_begin();
       it != end; ++it)
    {
      Elem *elem = *it;
      assert (elem);
      unsigned int elem_procid = elem->processor_id();

      // Assume we're partitioned before we try any AMR
      assert(elem_procid != DofObject::invalid_processor_id);

      ghost_elems_from_proc[elem_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_ids(libMesh::n_processors());

  // We know how many ghost elements live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      requested_ids[p].reserve(ghost_elems_from_proc[p]);

  for (MeshBase::element_iterator it = _mesh.elements_begin();
       it != end; ++it)
    {
      Elem *elem = *it;
      unsigned int elem_procid = elem->processor_id();

      requested_ids[elem_procid].push_back(elem->id());
    }

  // Set refinement flags from other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill;
      Parallel::send_receive(procup, requested_ids[procup],
                             procdown, request_to_fill);

      // Fill those requests
      std::vector<unsigned char> rflags(request_to_fill.size()),
                                 pflags(request_to_fill.size());
      for (unsigned int i=0; i != request_to_fill.size(); ++i)
        {
          Elem *elem = _mesh.elem(request_to_fill[i]);
          rflags[i] = elem->refinement_flag();
          pflags[i] = elem->p_refinement_flag();
        }

      // Trade back the results
      std::vector<unsigned char> ghost_rflags, ghost_pflags;
      Parallel::send_receive(procdown, rflags,
                             procup, ghost_rflags);
      Parallel::send_receive(procdown, pflags,
                             procup, ghost_pflags);
      assert (ghost_rflags.size() == requested_ids[procup].size());
      assert (ghost_pflags.size() == requested_ids[procup].size());

      // And see if we need to change any flags
      for (unsigned int i=0; i != requested_ids[procup].size(); ++i)
        {
          Elem *elem = _mesh.elem(requested_ids[procup][i]);
          unsigned char old_r_flag = elem->refinement_flag();
          unsigned char old_p_flag = elem->p_refinement_flag();
          if (old_r_flag != ghost_rflags[i])
            {
              elem->set_refinement_flag
                (static_cast<Elem::RefinementState>(ghost_rflags[i]));
              parallel_consistent = false;
            }
          if (old_p_flag != ghost_pflags[i])
            {
              elem->set_p_refinement_flag
                (static_cast<Elem::RefinementState>(ghost_pflags[i]));
              parallel_consistent = false;
            }
        }
    }

  // If we weren't consistent on any processor then we weren't
  // globally consistent
  Parallel::min(parallel_consistent);

  STOP_LOG ("make_flags_parallel_consistent()", "MeshRefinement");

  return parallel_consistent;
}



bool MeshRefinement::make_coarsening_compatible(const bool maintain_level_one)
{
  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("make_coarsening_compatible()", "MeshRefinement");
  
  bool _maintain_level_one = maintain_level_one;

  // If the user used non-default parameters, let's warn that they're
  // deprecated
  if (!maintain_level_one)
    {
      deprecated();
    }
  else
    _maintain_level_one = _face_level_mismatch_limit;

  
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
  
  // if there are no refined elements on this processor then
  // there is no work for us to do   
  if (max_level == 0 && max_p_level == 0)
    {
      STOP_LOG ("make_coarsening_compatible()", "MeshRefinement");

      // But we still have to check with other processors
      Parallel::min(compatible_with_refinement);

      return compatible_with_refinement;
    }

  
  
  // Loop over all the active elements.  If an element is marked
  // for coarsening we better check its neighbors.  If ANY of these neighbors
  // are marked for refinement AND are at the same level then there is a
  // conflict.  By convention refinement wins, so we un-mark the element for
  // coarsening.  Level-one would be violated in this case so we need to re-run
  // the loop.   
  if (_maintain_level_one)
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
              bool my_flag_changed = false;
	      
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
                              my_flag_changed = true;
                              break;
			    }
			}
		      else // I have a neighbor and it is not active. That means it has children.
			{  // While it _may_ be possible to coarsen us if all the children of
			   // that element want to be coarsened, it is impossible to know at this
			   // stage.  Forget about it for the moment...  This can be handled in
			   // two steps.
			  elem->set_refinement_flag(Elem::DO_NOTHING);
			  my_flag_changed = true;
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
			      my_flag_changed = true;
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
			             my_flag_changed = true;
                                     break;
			           }
                             }
                           if (my_flag_changed)
                             break;
			}
		}
		  
	      // If the current element's flag changed, we hadn't
	      // satisfied the level one rule.
	      if (my_flag_changed)
                level_one_satisfied = false;
	      
	      // Additionally, if it has non-local neighbors, and
	      // we're not in serial, then we'll eventually have to
	      // return compatible_with_refinement = false, because
	      // our change has to propagate to neighboring
	      // processors.
	      if (my_flag_changed && !_mesh.is_serial())
              for (unsigned int n=0; n != elem->n_neighbors(); ++n)
                {
                  Elem *neigh = elem->neighbor(n);
                  if (!neigh)
                    continue;
                  if (neigh->processor_id() !=
                      libMesh::processor_id())
                    {
                      compatible_with_refinement = false;
                      break;
                    }
                  // FIXME - for non-level one meshes we should
                  // test all descendants
                  if (neigh->has_children())
                    for (unsigned int c=0; c != neigh->n_children(); ++c)
                      if (neigh->child(c)->processor_id() !=
                          libMesh::processor_id())
                        {
                          compatible_with_refinement = false;
                          break;
                        }
                }
	    }      
	}
      while (!level_one_satisfied);
      
    } // end if (_maintain_level_one)
  
  
  // Next we look at all of the ancestor cells.
  // If there is a parent cell with all of its children
  // wanting to be unrefined then the element is a candidate
  // for unrefinement.  If all the children don't
  // all want to be unrefined then ALL of them need to have their
  // unrefinement flags cleared.    
  for (int level=(max_level); level >= 0; level--)
    {
      MeshBase::element_iterator       el     = _mesh.level_elements_begin(level);
      const MeshBase::element_iterator end_el = _mesh.level_elements_end(level); 
      for (; el != end_el; ++el)
        {
          Elem *elem = *el;
          if (elem->ancestor())
	    {
	  
	      // right now the element hasn't been disqualified
	      // as a candidate for unrefinement	   
	      bool is_a_candidate = true;
	      
	      for (unsigned int c=0; c<elem->n_children(); c++)
	        if ((elem->child(c)->refinement_flag() != Elem::COARSEN) ||
		    !elem->child(c)->active() )
	          is_a_candidate = false;
	      
	      if (!is_a_candidate)
	        {
	          elem->set_refinement_flag(Elem::INACTIVE);
		  
	          for (unsigned int c=0; c<elem->n_children(); c++)
		    if (elem->child(c)->refinement_flag() == Elem::COARSEN)
		      {
		        level_one_satisfied = false;
		        elem->child(c)->set_refinement_flag(Elem::DO_NOTHING);
		      }
	        }
	    }
         }
     }
     
  if (!level_one_satisfied && _maintain_level_one) goto repeat;
  
  
  // If all the children of a parent are set to be coarsened
  // then flag the parent so that they can kill thier kids...   
  MeshBase::element_iterator       all_el     = _mesh.elements_begin();
  const MeshBase::element_iterator all_el_end = _mesh.elements_end(); 
  for (; all_el != all_el_end; ++all_el)
    {
      Elem *elem = *all_el;
        if (elem->ancestor())
          {
	
	    // Presume all the children are flagged for coarsening
	    // and then look for a contradiction	 
	    bool all_children_flagged_for_coarsening = true;
	
	    for (unsigned int c=0; c<elem->n_children(); c++)
	      if (elem->child(c)->refinement_flag() != Elem::COARSEN)
	        all_children_flagged_for_coarsening = false;
	
	    if (all_children_flagged_for_coarsening)
	      elem->set_refinement_flag(Elem::COARSEN_INACTIVE);
            else
	      elem->set_refinement_flag(Elem::INACTIVE);
          }
    }
	
  STOP_LOG ("make_coarsening_compatible()", "MeshRefinement");

  // If one processor finds an incompatibility, we're globally
  // incompatible
  Parallel::min(compatible_with_refinement);
  
  return compatible_with_refinement;
}








bool MeshRefinement::make_refinement_compatible(const bool maintain_level_one)
{
  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("make_refinement_compatible()", "MeshRefinement");
  
  bool _maintain_level_one = maintain_level_one;

  // If the user used non-default parameters, let's warn that they're
  // deprecated
  if (!maintain_level_one)
    {
      deprecated();
    }
  else
    _maintain_level_one = _face_level_mismatch_limit;

  // Unless we encounter a specific situation level-one
  // will be satisfied after executing this loop just once
  bool level_one_satisfied = true;

  // Unless we encounter a specific situation we will be compatible
  // with any selected coarsening flags
  bool compatible_with_coarsening = true;

  // This loop enforces the level-1 rule.  We should only
  // execute it if the user indeed wants level-1 satisfied!   
  if (_maintain_level_one)
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
    } // end if (_maintain_level_one)

  // If we're not compatible on one processor, we're globally not
  // compatible
  Parallel::min(compatible_with_coarsening);
  
  STOP_LOG ("make_refinement_compatible()", "MeshRefinement");

  return compatible_with_coarsening;
}




bool MeshRefinement::_coarsen_elements ()
{
  // This function must be run on all processors at once
  parallel_only();

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

  // If the mesh changed on any processor, it changed globally
  Parallel::max(mesh_changed);
  
  STOP_LOG ("_coarsen_elements()", "MeshRefinement");

  return mesh_changed;
}







bool MeshRefinement::_refine_elements ()
{
  // This function must be run on all processors at once
  parallel_only();

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
  // To ensure that the new local nodes we add are given correct
  // processor ids, with ParallelMesh we'd better be adding elements
  // in increasing processor id order.  The default element sorting
  // should get that right.
#ifdef DEBUG
  unsigned int proc_id = 0;
#endif
  for (unsigned int e=0; e<local_copy_of_elements.size(); ++e)
    {
#ifdef DEBUG
      unsigned int next_proc_id =
        local_copy_of_elements[e]->processor_id();
      assert (_mesh.is_serial() || next_proc_id >= proc_id);
      proc_id = next_proc_id;
#endif
      local_copy_of_elements[e]->refine(*this);
    }

  if (!_mesh.is_serial())
    {
      this->make_nodes_parallel_consistent();
      this->make_elems_parallel_consistent();
    }
  
  // Clear the _new_nodes_map and _unused_elements data structures.
  this->clear();
  
  // If the mesh changed on any processor, it changed globally
  Parallel::max(mesh_changed);
  
  STOP_LOG ("_refine_elements()", "MeshRefinement");

  return mesh_changed;
}



void MeshRefinement::make_nodes_parallel_consistent()
{
  // This function must be run on all processors at once
  parallel_only();

  // Local nodes in the _new_nodes_map have authoritative ids
  // and correct processor ids, nodes touching local elements
  // have correct processor ids but need to have their ids
  // corrected, and ghost nodes not touching local elements
  // may need to have both ids and processor ids corrected.

  // First correct the processor ids, so we'll ask the 
  // right processor when correcting ids later

  // Count the nodes to ask each processor about
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  const MeshBase::node_iterator end = _mesh.nodes_end();

  for (MeshBase::node_iterator it  = _mesh.nodes_begin();
       it != end; ++it)
    {
      Node *node = *it;
      assert (node);
      unsigned int node_procid = node->processor_id();
      assert (node_procid != DofObject::invalid_processor_id);

      ghost_objects_from_proc[node_procid]++;
    }

  // Request sets to send to each processor on the first pass
  std::vector<std::vector<Real> >
    requested_nodes_x(libMesh::n_processors()),
    requested_nodes_y(libMesh::n_processors()),
    requested_nodes_z(libMesh::n_processors());
  // Corresponding temporary ids to keep track of
  std::vector<std::vector<unsigned int> >
    requested_nodes_id(libMesh::n_processors());
  // And (hopefully much smaller) sets for the second pass
  std::vector<std::vector<Real> >
    rerequested_nodes_x(libMesh::n_processors()),
    rerequested_nodes_y(libMesh::n_processors()),
    rerequested_nodes_z(libMesh::n_processors());
  std::vector<std::vector<unsigned int> >
    rerequested_nodes_id(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        requested_nodes_x[p].reserve(ghost_objects_from_proc[p]);
        requested_nodes_y[p].reserve(ghost_objects_from_proc[p]);
        requested_nodes_z[p].reserve(ghost_objects_from_proc[p]);
        requested_nodes_id[p].reserve(ghost_objects_from_proc[p]);
      }

  for (MeshBase::node_iterator it  = _mesh.nodes_begin();
       it != end; ++it)
    {
      Node *node = *it;
      unsigned int node_procid = node->processor_id();

      requested_nodes_x[node_procid].push_back((*node)(0));
      requested_nodes_y[node_procid].push_back((*node)(1));
      requested_nodes_z[node_procid].push_back((*node)(2));
      requested_nodes_id[node_procid].push_back(node->id());
    }

  // Trade requests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<Real> request_to_fill_x,
                        request_to_fill_y,
                        request_to_fill_z;
      Parallel::send_receive(procup, requested_nodes_x[procup],
                             procdown, request_to_fill_x);
      Parallel::send_receive(procup, requested_nodes_y[procup],
                             procdown, request_to_fill_y);
      Parallel::send_receive(procup, requested_nodes_z[procup],
                             procdown, request_to_fill_z);
      assert (request_to_fill_x.size() == request_to_fill_y.size());
      assert (request_to_fill_x.size() == request_to_fill_z.size());

      // Find the processor id (and if it's local, the id)
      // of each requested node
      std::vector<unsigned int> node_proc_ids(request_to_fill_x.size()),
                                node_ids(request_to_fill_x.size());
      for (unsigned int i=0; i != request_to_fill_x.size(); ++i)
        {
          Point p(request_to_fill_x[i],
                  request_to_fill_y[i],
                  request_to_fill_z[i]);
          unsigned int key = this->point_key(p);

          // Look for this point in the multimap
          std::pair<map_type::iterator, map_type::iterator>
            pos = _new_nodes_map.equal_range(key);

          // We'd better find every node we're asked for
          assert (pos.first != pos.second);

          Node *node = NULL;
          // FIXME - what tolerance should we use?
          for (; pos.first != pos.second; ++pos.first)
            if (p.absolute_fuzzy_equals(*(pos.first->second), TOLERANCE))
              {
                node = pos.first->second;
                break;
              }
            else
              {
                // Make sure this map conflict isn't a key bug
                assert (this->point_key(*(pos.first->second)) == 
                        this->point_key(p));
              }

          // We'd better have found every node we're asked for
          assert (node);

          // Return the node's correct processor id,
          // and our (correct if it's local) id for it.
          node_proc_ids[i] = node->processor_id();
          node_ids[i] = node->id();
        }
      
      // Trade back the results
      std::vector<unsigned int> filled_node_proc_ids, filled_node_ids;
      Parallel::send_receive(procdown, node_proc_ids,
                             procup, filled_node_proc_ids);
      Parallel::send_receive(procdown, node_ids,
                             procup, filled_node_ids);
      assert (requested_nodes_x[procup].size() == filled_node_proc_ids.size());
      assert (requested_nodes_x[procup].size() == filled_node_ids.size());

      // Set the ghost node processor ids and ids we've now been
      // informed of, and build request sets for ids we need to
      // request from a different processor
      for (unsigned int i=0; i != filled_node_ids.size(); ++i)
        {
          Node *node = _mesh.node_ptr(requested_nodes_id[procup][i]);
          const unsigned int new_procid = filled_node_proc_ids[i];

          // Set ids of and move nodes where we asked their local processor
          assert (node->processor_id() == procup);
          if (procup == new_procid)
            {
              const unsigned int old_id = requested_nodes_id[procup][i],
                                 new_id = filled_node_ids[i];
              if (old_id != new_id)
                _mesh.renumber_node(old_id, new_id);
            }

          // Rerequest ids of nodes where we should have asked another
          // processor
          else
            {
              node->processor_id() = new_procid;
              // We need to rerequest this node's id, from the
              // correct processor this time.

              // There's no obvious way to reserve() this vector,
              // so let's hope our STL amortizes resize() properly.
              // O(N) overhead on small vectors shouldn't be bad.
              rerequested_nodes_x[new_procid].push_back((*node)(0));
              rerequested_nodes_y[new_procid].push_back((*node)(1));
              rerequested_nodes_z[new_procid].push_back((*node)(2));
              rerequested_nodes_id[new_procid].push_back(node->id());
            }
        }
    }

  // We're done with the first set of requests
  requested_nodes_x.clear();
  requested_nodes_y.clear();
  requested_nodes_z.clear();
  requested_nodes_id.clear();

  // Trade rerequests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my rerequests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<Real> rerequest_to_fill_x,
                        rerequest_to_fill_y,
                        rerequest_to_fill_z;
      Parallel::send_receive(procup, rerequested_nodes_x[procup],
                             procdown, rerequest_to_fill_x);
      Parallel::send_receive(procup, rerequested_nodes_y[procup],
                             procdown, rerequest_to_fill_y);
      Parallel::send_receive(procup, rerequested_nodes_z[procup],
                             procdown, rerequest_to_fill_z);
      assert (rerequest_to_fill_x.size() == rerequest_to_fill_y.size());
      assert (rerequest_to_fill_x.size() == rerequest_to_fill_z.size());

      // Find the id of each rerequested node
      std::vector<unsigned int> node_ids(rerequest_to_fill_x.size());
      for (unsigned int i=0; i != rerequest_to_fill_x.size(); ++i)
        {
          Point p(rerequest_to_fill_x[i],
                  rerequest_to_fill_y[i],
                  rerequest_to_fill_z[i]);
          unsigned int key = this->point_key(p);

          // Look for this point in the multimap
          std::pair<map_type::iterator, map_type::iterator>
            pos = _new_nodes_map.equal_range(key);

          // We'd better find every node we're asked for
          assert (pos.first != pos.second);

          Node *node = NULL;
          // FIXME - what tolerance should we use?
          for (; pos.first != pos.second; ++pos.first)
            if (p.absolute_fuzzy_equals(*(pos.first->second), TOLERANCE))
              {
                node = pos.first->second;
                break;
              }

          // We'd better have found every node we're asked for
          assert (node);

          // Return the node's correct id
          node_ids[i] = node->id();
        }
      
      // Trade back the results
      std::vector<unsigned int> filled_node_ids;
      Parallel::send_receive(procdown, node_ids,
                             procup, filled_node_ids);
      assert (rerequested_nodes_x[procup].size() == filled_node_ids.size());

      // Set the ghost node ids we've now been informed of
      for (unsigned int i=0; i != filled_node_ids.size(); ++i)
        {
          const unsigned int old_id = rerequested_nodes_id[procup][i],
                             new_id = filled_node_ids[i];
          if (old_id != new_id)
            _mesh.renumber_node(old_id, new_id);
        }
    }
}


void MeshRefinement::make_elems_parallel_consistent()
{
  // This function must be run on all processors at once
  parallel_only();

  // Newly added local elements have authoritative ids
  // and correct processor ids, ghost elements have correct
  // processor ids but need to have their ids corrected.

  // Count the elements to ask each processor about
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  const MeshBase::element_iterator end = _mesh.elements_end();

  for (MeshBase::element_iterator it  = _mesh.elements_begin();
       it != end; ++it)
    {
      Elem *elem = *it;
      unsigned int elem_procid = elem->processor_id();

      // We only have to worry about new ghost child elements
      if (!elem->parent() || !elem->active() ||
          elem_procid == libMesh::processor_id())
        continue;

      assert (elem_procid != DofObject::invalid_processor_id);

      ghost_objects_from_proc[elem_procid]++;
    }

  // Request sets to send to each processor on the first pass
  std::vector<std::vector<unsigned int> >
    requested_parent_ids(libMesh::n_processors()),
    requested_child_nums(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        requested_parent_ids[p].reserve(ghost_objects_from_proc[p]);
        requested_child_nums[p].reserve(ghost_objects_from_proc[p]);
      }

  for (MeshBase::element_iterator it  = _mesh.elements_begin();
       it != end; ++it)
    {
      Elem *elem = *it;
      unsigned int elem_procid = elem->processor_id();
      const Elem *parent = elem->parent();

      // We only have to worry about new ghost child elements
      if (!parent || !elem->active() ||
          elem_procid == libMesh::processor_id())
        continue;

      requested_parent_ids[elem_procid].push_back(parent->id());
      requested_child_nums[elem_procid].push_back
        (parent->which_child_am_i(elem));
    }
#ifdef DEBUG
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        assert(requested_parent_ids[p].size() == ghost_objects_from_proc[p]);
        assert(requested_child_nums[p].size() == ghost_objects_from_proc[p]);
      }
#endif

  // Trade requests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill_parent_ids,
                                request_to_fill_child_nums;
      Parallel::send_receive(procup, requested_parent_ids[procup],
                             procdown, request_to_fill_parent_ids);
      Parallel::send_receive(procup, requested_child_nums[procup],
                             procdown, request_to_fill_child_nums);
      assert (request_to_fill_parent_ids.size() ==
              request_to_fill_child_nums.size());

      // Find the id of each requested element
      std::vector<unsigned int> elem_ids(request_to_fill_parent_ids.size());
      for (unsigned int i=0; i != request_to_fill_parent_ids.size(); ++i)
        {
          Elem *parent = _mesh.elem(request_to_fill_parent_ids[i]);
          assert (parent);
          assert (parent->has_children());
          Elem *child = parent->child(request_to_fill_child_nums[i]);
          assert (child);
          assert (child->active());

          // Return the child element's correct id
          elem_ids[i] = child->id();
        }
      
      // Trade back the results
      std::vector<unsigned int> filled_elem_ids;
      Parallel::send_receive(procdown, elem_ids,
                             procup, filled_elem_ids);
      assert (requested_parent_ids[procup].size() == filled_elem_ids.size());

      // Set those ghost element ids
      for (unsigned int i=0; i != filled_elem_ids.size(); ++i)
        {
          Elem *parent = _mesh.elem(requested_parent_ids[procup][i]);
          assert (parent);
          assert (parent->has_children());
          Elem *child = parent->child(requested_child_nums[procup][i]);
          assert (child);
          assert (child->active());
          const unsigned int old_id = child->id(),
                             new_id = filled_elem_ids[i];
          if (old_id != new_id)
            _mesh.renumber_elem(old_id, new_id);
        }
    }
}



unsigned int MeshRefinement::point_key (const Point &p) const
{
  Real xscaled = (p(0) - _lower_bound[0])/
                 (_upper_bound[0] - _lower_bound[0]),
       yscaled = (p(1) - _lower_bound[1])/
                 (_upper_bound[1] - _lower_bound[1]),
       zscaled = (p(2) - _lower_bound[2])/
                 (_upper_bound[2] - _lower_bound[2]);

  // 10 bits per coordinate, to work with 32+ bit machines
  unsigned chunkmax = 1024;
  Real chunkfloat = 1024.0;

  unsigned int n0 = static_cast<unsigned int> (chunkfloat * xscaled),
               n1 = static_cast<unsigned int> (chunkfloat * yscaled),
               n2 = static_cast<unsigned int> (chunkfloat * zscaled);

  return chunkmax*chunkmax*n0 + chunkmax*n1 + n2;
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
