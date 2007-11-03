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
#include "partitioner.h"
#include "mesh_base.h"
#include "elem.h"


// ------------------------------------------------------------
// Partitioner implementation


void Partitioner::partition (MeshBase& mesh,
			     const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;

  // Call the partitioning function
  this->_do_partition(mesh,n);

  // Set the node's processor ids
  this->_set_node_processor_ids(mesh);
}





void Partitioner::repartition (MeshBase& mesh,
			       const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;
  
  // Call the partitioning function
  this->_do_repartition(mesh,n);

  // Set the node's processor ids
  this->_set_node_processor_ids(mesh);
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




void Partitioner::_set_node_processor_ids(MeshBase& mesh)
{
  // Unset any previously-set node processor ids
  // (maybe from previous partitionings).
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();
  
  for ( ; node_it != node_end; ++node_it)
    (*node_it)->invalidate_processor_id();
  
  
  // Loop over all the elements
  MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 
  
  for ( ; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;
      
      // For each node, set the processor ID to the min of
      // its current value and this Element's processor id.
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
	elem->get_node(n)->processor_id() = std::min(elem->get_node(n)->processor_id(),
						     elem->processor_id());
    }
}
