// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "elem.h"
#include "mesh_base.h"
#include "parallel.h"
#include "partitioner.h"
#include "mesh_tools.h"
#include "mesh_communication.h"
#include "libmesh_logging.h"

namespace libMesh
{



// ------------------------------------------------------------
// Partitioner static data
const unsigned int Partitioner::communication_blocksize = 1000000;



// ------------------------------------------------------------
// Partitioner implementation
void Partitioner::partition (MeshBase& mesh,
			     const unsigned int n)
{
  parallel_only();

  // BSK - temporary fix while redistribution is integrated 6/26/2008
  // Uncomment this to not repartition in parallel
//   if (!mesh.is_serial())
//     return;

  // we cannot partition into more pieces than we have
  // active elements!
  const unsigned int n_parts =
    std::min(mesh.n_active_elem(), n);
  
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n_parts;

  if (n_parts == 1)
    {
      this->single_partition (mesh);
      return;
    }
  
  // First assign a temporary partitioning to any unpartitioned elements
  Partitioner::partition_unpartitioned_elements(mesh, n_parts);
  
  // Call the partitioning function
  this->_do_partition(mesh,n_parts);

  // Set the parent's processor ids
  Partitioner::set_parent_processor_ids(mesh);

  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);

#ifdef DEBUG
  // Messed up elem processor_id()s can leave us without the child
  // elements we need to restrict vectors on a distributed mesh
  MeshTools::libmesh_assert_valid_elem_procids(mesh);

  MeshTools::libmesh_assert_valid_node_procids(mesh);
#endif
}





void Partitioner::repartition (MeshBase& mesh,
			       const unsigned int n)
{
  // we cannot partition into more pieces than we have
  // active elements!
  const unsigned int n_parts =
    std::min(mesh.n_active_elem(), n);
  
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n_parts;

  if (n_parts == 1)
    {
      this->single_partition (mesh);
      return;
    }
  
  // First assign a temporary partitioning to any unpartitioned elements
  Partitioner::partition_unpartitioned_elements(mesh, n_parts);
  
  // Call the partitioning function
  this->_do_repartition(mesh,n_parts);
  
  // Set the parent's processor ids
  Partitioner::set_parent_processor_ids(mesh);
  
  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);
}





void Partitioner::single_partition (MeshBase& mesh)
{
  START_LOG("single_partition()","Partitioner");
  
  // Loop over all the elements and assign them to processor 0.
  MeshBase::element_iterator       elem_it  = mesh.elements_begin();
  const MeshBase::element_iterator elem_end = mesh.elements_end(); 

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->processor_id() = 0;

  // For a single partition, all the nodes are on processor 0
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();
  
  for ( ; node_it != node_end; ++node_it)
    (*node_it)->processor_id() = 0;

  STOP_LOG("single_partition()","Partitioner");
}



void Partitioner::partition_unpartitioned_elements (MeshBase &mesh,
						    const unsigned int n_subdomains)
{
  MeshBase::const_element_iterator       it  = mesh.unpartitioned_elements_begin();
  const MeshBase::const_element_iterator end = mesh.unpartitioned_elements_end();

  const unsigned int n_unpartitioned_elements = MeshTools::n_elem (it, end);

  // the unpartitioned elements must exist on all processors. If the range is empty on one
  // it is empty on all, and we can quit right here.
  if (!n_unpartitioned_elements) return;
	 
  // find the target subdomain sizes
  std::vector<unsigned int> subdomain_bounds(libMesh::n_processors());

  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
    {
      unsigned int tgt_subdomain_size = 0;

      // watch out for the case that n_subdomains < n_processors
      if (pid < n_subdomains)
	{
	  tgt_subdomain_size = n_unpartitioned_elements/n_subdomains;
      
	  if (pid < n_unpartitioned_elements%n_subdomains)
	    tgt_subdomain_size++;

	}
      
      //libMesh::out << "pid, #= " << pid << ", " << tgt_subdomain_size << std::endl;
      if (pid == 0)
	subdomain_bounds[0] = tgt_subdomain_size;
      else
	subdomain_bounds[pid] = subdomain_bounds[pid-1] + tgt_subdomain_size;
    }
  
  libmesh_assert (subdomain_bounds.back() == n_unpartitioned_elements);  
  
  // create the unique mapping for all unpartitioned elements independent of partitioning
  // determine the global indexing for all the unpartitoned elements
  std::vector<unsigned int> global_indices;
    
  // Calling this on all processors a unique range in [0,n_unpartitioned_elements) is constructed.  
  // Only the indices for the elements we pass in are returned in the array.
  MeshCommunication().find_global_indices (MeshTools::bounding_box(mesh), it, end, 
					   global_indices);
  
  for (unsigned int cnt=0; it != end; ++it)
    {
      Elem *elem = *it;
      
      libmesh_assert (cnt < global_indices.size());
      const unsigned int global_index =
	global_indices[cnt++];
      
      libmesh_assert (global_index < subdomain_bounds.back());
      libmesh_assert (global_index < n_unpartitioned_elements);

      const unsigned int subdomain_id =
	std::distance(subdomain_bounds.begin(),
		      std::upper_bound(subdomain_bounds.begin(),
				       subdomain_bounds.end(),
				       global_index));
      libmesh_assert (subdomain_id < n_subdomains);
     
      elem->processor_id() = subdomain_id;		
      //libMesh::out << "assigning " << global_index << " to " << subdomain_id << std::endl;	      
    }
}



void Partitioner::set_parent_processor_ids(MeshBase& mesh)
{
  START_LOG("set_parent_processor_ids()","Partitioner");
  
#ifdef LIBMESH_ENABLE_AMR

  // If the mesh is serial we have access to all the elements,
  // in particular all the active ones.  We can therefore set
  // the parent processor ids indirecly through their children, and
  // set the subactive processor ids while examining their active
  // ancestors.
  // By convention a parent is assigned to the minimum processor
  // of all its children, and a subactive is assigned to the processor
  // of its active ancestor.
  if (mesh.is_serial())
    {
      // Loop over all the active elements in the mesh
      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end();
      
      for ( ; it!=end; ++it)
	{
	  Elem *child  = *it;

          // First set descendents

          std::vector<const Elem*> subactive_family;
          child->total_family_tree(subactive_family);
          for (unsigned int i = 0; i != subactive_family.size(); ++i)
            const_cast<Elem*>(subactive_family[i])->processor_id() = child->processor_id();

          // Then set ancestors

	  Elem *parent = child->parent();

	  while (parent)
	    {
	      // invalidate the parent id, otherwise the min below
	      // will not work if the current parent id is less
	      // than all the children!
	      parent->invalidate_processor_id();
	      
	      for(unsigned int c=0; c<parent->n_children(); c++)
		{
		  child = parent->child(c);
		  libmesh_assert(child);
		  libmesh_assert(!child->is_remote());
		  libmesh_assert(child->processor_id() != DofObject::invalid_processor_id);
		  parent->processor_id() = std::min(parent->processor_id(),
						    child->processor_id());
		}	      
	      parent = parent->parent();
	    }
	}
    }

  // When the mesh is parallel we cannot guarantee that parents have access to
  // all their children.
  else
    {
      // Setting subactive processor ids is easy: we can guarantee
      // that children have access to all their parents.

      // Loop over all the active elements in the mesh
      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end();
      
      for ( ; it!=end; ++it)
        {
	  Elem *child  = *it;

          std::vector<const Elem*> subactive_family;
          child->total_family_tree(subactive_family);
          for (unsigned int i = 0; i != subactive_family.size(); ++i)
            const_cast<Elem*>(subactive_family[i])->processor_id() = child->processor_id();
        }

      // When the mesh is parallel we cannot guarantee that parents have access to
      // all their children.

      // We will use a brute-force approach here.  Each processor finds its parent
      // elements and sets the parent pid to the minimum of its local children.
      // A global reduction is then performed to make sure the true minimum is found.
      // As noted, this is required because we cannot guarantee that a parent has
      // access to all its children on any single processor.
      parallel_only();
      libmesh_assert(MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
			       mesh.unpartitioned_elements_end()) == 0);

      const unsigned int max_elem_id = mesh.max_elem_id();

      std::vector<processor_id_type>
	parent_processor_ids (std::min(communication_blocksize,
				       max_elem_id));
      
      for (unsigned int blk=0, last_elem_id=0; last_elem_id<max_elem_id; blk++)
	{
 	                      last_elem_id = std::min((blk+1)*communication_blocksize, max_elem_id);
	  const unsigned int first_elem_id = blk*communication_blocksize;

	  std::fill (parent_processor_ids.begin(),
		     parent_processor_ids.end(),
		     DofObject::invalid_processor_id);

	  // first build up local contributions to parent_processor_ids
	  MeshBase::element_iterator       not_it  = mesh.ancestor_elements_begin();
	  const MeshBase::element_iterator not_end = mesh.ancestor_elements_end(); 

	  bool have_parent_in_block = false;
	  
	  for ( ; not_it != not_end; ++not_it)
	    {
	      Elem *parent = *not_it;

	      const unsigned int parent_idx = parent->id();
	      libmesh_assert (parent_idx < max_elem_id);

	      if ((parent_idx >= first_elem_id) &&
		  (parent_idx <  last_elem_id))
		{
		  have_parent_in_block = true;
		  processor_id_type parent_pid = DofObject::invalid_processor_id;

		  for (unsigned int c=0; c<parent->n_children(); c++)
		    parent_pid = std::min (parent_pid, parent->child(c)->processor_id());
		  
		  const unsigned int packed_idx = parent_idx - first_elem_id;
		  libmesh_assert (packed_idx < parent_processor_ids.size());

		  parent_processor_ids[packed_idx] = parent_pid;
		}
	    }

	  // then find the global minimum
	  Parallel::min (parent_processor_ids);

	  // and assign the ids, if we have a parent in this block.
	  if (have_parent_in_block)
	    for (not_it = mesh.ancestor_elements_begin();
		 not_it != not_end; ++not_it)
	      {
		Elem *parent = *not_it;
		
		const unsigned int parent_idx = parent->id();
		
		if ((parent_idx >= first_elem_id) &&
		    (parent_idx <  last_elem_id))
		  {
		    const unsigned int packed_idx = parent_idx - first_elem_id;
		    libmesh_assert (packed_idx < parent_processor_ids.size());
		    
		    const processor_id_type parent_pid =
		      parent_processor_ids[packed_idx];
		    
		    libmesh_assert (parent_pid != DofObject::invalid_processor_id);
		    
		    parent->processor_id() = parent_pid;
		  }
	      }
	}
    }
  
#endif // LIBMESH_ENABLE_AMR

  STOP_LOG("set_parent_processor_ids()","Partitioner");
}



void Partitioner::set_node_processor_ids(MeshBase& mesh)
{
  START_LOG("set_node_processor_ids()","Partitioner");

  // This function must be run on all processors at once
  parallel_only();

  // If we have any unpartitioned elements at this 
  // stage there is a problem
  libmesh_assert (MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
			    mesh.unpartitioned_elements_end()) == 0);


//   const unsigned int orig_n_local_nodes = mesh.n_local_nodes();

//   libMesh::err << "[" << libMesh::processor_id() << "]: orig_n_local_nodes="
// 	    << orig_n_local_nodes << std::endl;

  // Build up request sets.  Each node is currently owned by a processor because
  // it is connected to an element owned by that processor.  However, during the
  // repartitioning phase that element may have been assigned a new processor id, but
  // it is still resident on the original processor.  We need to know where to look
  // for new ids before assigning new ids, otherwise we may be asking the wrong processors
  // for the wrong information.
  //
  // The only remaining issue is what to do with unpartitioned nodes.  Since they are required
  // to live on all processors we can simply rely on ourselves to number them properly.
  std::vector<std::vector<unsigned int> >
    requested_node_ids(libMesh::n_processors());

  // Loop over all the nodes, count the ones on each processor.  We can skip ourself
  std::vector<unsigned int> ghost_nodes_from_proc(libMesh::n_processors(), 0);

  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();
  
  for (; node_it != node_end; ++node_it)
    {
      Node *node = *node_it;
      libmesh_assert(node);
      const unsigned int current_pid = node->processor_id();
      if (current_pid != libMesh::processor_id() &&
	  current_pid != DofObject::invalid_processor_id)
	{
	  libmesh_assert(current_pid < ghost_nodes_from_proc.size());
	  ghost_nodes_from_proc[current_pid]++;
	}
    }

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int pid=0; pid != libMesh::n_processors(); ++pid)
    requested_node_ids[pid].reserve(ghost_nodes_from_proc[pid]);

  // We need to get the new pid for each node from the processor
  // which *currently* owns the node.  We can safely skip ourself
  for (node_it = mesh.nodes_begin(); node_it != node_end; ++node_it)
    {
      Node *node = *node_it;
      libmesh_assert(node);
      const processor_id_type current_pid = node->processor_id();      
      if (current_pid != libMesh::processor_id() &&
	  current_pid != DofObject::invalid_processor_id)
	{
	  libmesh_assert(current_pid < requested_node_ids.size());
	  libmesh_assert(requested_node_ids[current_pid].size() <
		 ghost_nodes_from_proc[current_pid]);
	  requested_node_ids[current_pid].push_back(node->id());
	}
      
      // Unset any previously-set node processor ids
      node->invalidate_processor_id();
    }
  
  // Loop over all the active elements
  MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 
  
  for ( ; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;
      libmesh_assert(elem);

      libmesh_assert (elem->processor_id() != DofObject::invalid_processor_id);
      
      // For each node, set the processor ID to the min of
      // its current value and this Element's processor id.
      // 
      // TODO: we would probably get better parallel partitioning if
      // we did something like "min for even numbered nodes, max for
      // odd numbered".  We'd need to be careful about how that would
      // affect solution ordering for I/O, though.
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
	elem->get_node(n)->processor_id() = std::min(elem->get_node(n)->processor_id(),
						     elem->processor_id());
    }

  // And loop over the subactive elements, but don't reassign
  // nodes that are already active on another processor.
  MeshBase::element_iterator       sub_it  = mesh.subactive_elements_begin();
  const MeshBase::element_iterator sub_end = mesh.subactive_elements_end(); 
  
  for ( ; sub_it != sub_end; ++sub_it)
    {
      Elem* elem = *sub_it;
      libmesh_assert(elem);

      libmesh_assert (elem->processor_id() != DofObject::invalid_processor_id);
      
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        if (elem->get_node(n)->processor_id() == DofObject::invalid_processor_id)
	  elem->get_node(n)->processor_id() = elem->processor_id();
    }

  // Same for the inactive elements -- we will have already gotten most of these
  // nodes, *except* for the case of a parent with a subset of children which are
  // ghost elements.  In that case some of the parent nodes will not have been
  // properly handled yet
  MeshBase::element_iterator       not_it  = mesh.not_active_elements_begin();
  const MeshBase::element_iterator not_end = mesh.not_active_elements_end(); 
  
  for ( ; not_it != not_end; ++not_it)
    {
      Elem* elem = *not_it;
      libmesh_assert(elem);

      libmesh_assert (elem->processor_id() != DofObject::invalid_processor_id);
      
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        if (elem->get_node(n)->processor_id() == DofObject::invalid_processor_id)
	  elem->get_node(n)->processor_id() = elem->processor_id();
    }

#ifndef NDEBUG
  {
    // make sure all the nodes connected to any element have received a
    // valid processor id
    std::set<const Node*> used_nodes;
    MeshBase::element_iterator       all_it  = mesh.elements_begin();
    const MeshBase::element_iterator all_end = mesh.elements_end(); 
  
    for ( ; all_it != all_end; ++all_it)
      {
	Elem* elem = *all_it;
	libmesh_assert(elem);
	libmesh_assert(elem->processor_id() != DofObject::invalid_processor_id);
	for (unsigned int n=0; n<elem->n_nodes(); ++n)
	  used_nodes.insert(elem->get_node(n));
      }

    for (node_it = mesh.nodes_begin(); node_it != node_end; ++node_it)
      {
	Node *node = *node_it;
	libmesh_assert(node);
	libmesh_assert(used_nodes.count(node));
	libmesh_assert(node->processor_id() != DofObject::invalid_processor_id);
      }
  }
#endif

  // Next set node ids from other processors, excluding self
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill;
      Parallel::send_receive(procup, requested_node_ids[procup],
                             procdown, request_to_fill);

      // Fill those requests in-place
      for (unsigned int i=0; i != request_to_fill.size(); ++i)
        {
          Node *node = mesh.node_ptr(request_to_fill[i]);
          libmesh_assert(node);
	  const processor_id_type new_pid = node->processor_id();
	  libmesh_assert (new_pid != DofObject::invalid_processor_id);
	  libmesh_assert (new_pid < mesh.n_partitions()); // this is the correct test --
          request_to_fill[i] = new_pid;           //  the number of partitions may
        }                                         //  not equal the number of processors

      // Trade back the results
      std::vector<unsigned int> filled_request;
      Parallel::send_receive(procdown, request_to_fill,
                             procup,   filled_request);
      libmesh_assert(filled_request.size() == requested_node_ids[procup].size());
      
      // And copy the id changes we've now been informed of
      for (unsigned int i=0; i != filled_request.size(); ++i)
        {
          Node *node = mesh.node_ptr(requested_node_ids[procup][i]);
	  libmesh_assert(node);
          libmesh_assert(filled_request[i] < mesh.n_partitions()); // this is the correct test --
          node->processor_id(filled_request[i]);           //  the number of partitions may
        }                                                  //  not equal the number of processors
    }

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_node_procids(mesh);
#endif
  
  STOP_LOG("set_node_processor_ids()","Partitioner");
}

} // namespace libMesh
