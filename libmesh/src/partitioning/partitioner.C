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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "elem.h"
#include "mesh_base.h"
#include "parallel.h"
#include "partitioner.h"


// ------------------------------------------------------------
// Partitioner implementation


void Partitioner::partition (MeshBase& mesh,
			     const unsigned int n)
{
  // For now we don't repartition in parallel
  if (!mesh.is_serial())
    return;

  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;

  // Call the partitioning function
  this->_do_partition(mesh,n);

  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);
}





void Partitioner::repartition (MeshBase& mesh,
			       const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;
  
  // Call the partitioning function
  this->_do_repartition(mesh,n);

  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);
}





void Partitioner::single_partition (MeshBase& mesh)
{
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
}




void Partitioner::set_node_processor_ids(MeshBase& mesh)
{
  // This function must be run on all processors at once
  parallel_only();

  // Unset any previously-set node processor ids
  // (maybe from previous partitionings).
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();
  
  for ( ; node_it != node_end; ++node_it)
    {
      Node *node = *node_it;
      assert(node);
      node->invalidate_processor_id();
    }
  
  
  // Loop over all the active elements
  MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 
  
  for ( ; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;
      assert(elem);
      
      // For each node, set the processor ID to the min of
      // its current value and this Element's processor id.
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
      assert(elem);
      
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        if (elem->get_node(n)->processor_id() == DofObject::invalid_processor_id)
	  elem->get_node(n)->processor_id() = std::min(elem->get_node(n)->processor_id(),
						       elem->processor_id());
    }

  // At this point, if we're in parallel, all our own processor ids
  // should be correct, but ghost nodes may be incorrect.  However,
  // our tenative processor ids will let us know who to ask.

  // Loop over all the nodes, count the ones on each processor
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  node_it  = mesh.nodes_begin();
  for ( ; node_it != node_end; ++node_it)
    {
      Node *node = *node_it;
      assert(node);
      unsigned int obj_procid = node->processor_id();
      assert(obj_procid != DofObject::invalid_processor_id);
      ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_ids(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      requested_ids[p].reserve(ghost_objects_from_proc[p]);

  node_it  = mesh.nodes_begin();
  for ( ; node_it != node_end; ++node_it)
    {
      Node *node = *node_it;
      requested_ids[node->processor_id()].push_back(node->id());
    }

  // Next set ghost object ids from other processors
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
      std::vector<unsigned int> new_ids(request_to_fill.size());
      for (unsigned int i=0; i != request_to_fill.size(); ++i)
        {
          Node *node = mesh.node_ptr(request_to_fill[i]);
          assert(node);
          new_ids[i] = node->processor_id();
        }

      // Trade back the results
      std::vector<unsigned int> filled_request;
      Parallel::send_receive(procdown, new_ids,
                             procup, filled_request);
      assert(filled_request.size() == requested_ids[procup].size());
      
      // And copy the id changes we've now been informed of
      for (unsigned int i=0; i != filled_request.size(); ++i)
        {
          Node *node = mesh.node_ptr(requested_ids[procup][i]);
          assert(filled_request[i] < libMesh::n_processors());
          node->processor_id(filled_request[i]);
        }
    }
}
