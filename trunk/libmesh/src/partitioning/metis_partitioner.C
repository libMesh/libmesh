// $Id: metis_partitioner.C,v 1.17 2004-11-14 18:51:59 jwpeterson Exp $

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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_config.h"
#include "mesh_base.h"
#include "metis_partitioner.h"
#include "libmesh_logging.h"
#include "elem.h"

#ifdef HAVE_METIS
  namespace Metis {
    extern "C" {
#     include "metis.h"
    }
  }
#else
#  include "sfc_partitioner.h"
#endif



// ------------------------------------------------------------
// MetisPartitioner implementation
void MetisPartitioner::_do_partition (MeshBase& mesh,
				      const unsigned int n_pieces)
{
  assert (n_pieces > 0);

  // Check for an easy return
  if (n_pieces == 1)
    {
      this->single_partition (mesh);
      return;
    }

// What to do if the Metis library IS NOT present
#ifndef HAVE_METIS

  here();
  std::cerr << "ERROR: The library has been built without"    << std::endl
	    << "Metis support.  Using a space-filling curve"  << std::endl
	    << "partitioner instead!"                         << std::endl;

  SFCPartitioner sfcp;

  sfcp.partition (mesh, n_pieces);
  
// What to do if the Metis library IS present
#else

  START_LOG("partition()", "MetisPartitioner");

  const unsigned int n_active_elem = mesh.n_active_elem();
  const unsigned int n_elem        = mesh.n_elem();
  
  // build the graph
  // the forward_map maps each active element id
  // into a contiguous block of indices for Metis
  std::vector<unsigned int> forward_map (n_elem, libMesh::invalid_uint);
  
  std::vector<int> xadj;
  std::vector<int> adjncy;
  std::vector<int> options(5);
  std::vector<int> vwgt(n_active_elem);
  std::vector<int> part(n_active_elem);

  xadj.reserve(n_active_elem+1);
  
  int
    n = static_cast<int>(n_active_elem),  // number of "nodes" (elements)
                                          //   in the graph
    wgtflag = 2,                          // weights on vertices only,
                                          //   none on edges
    numflag = 0,                          // C-style 0-based numbering
    nparts  = static_cast<int>(n_pieces), // number of subdomains to create
    edgecut = 0;                          // the numbers of edges cut by the
                                          //   resulting partition

  // Set the options
  options[0] = 0; // use default options


  // Metis will only consider the active elements.
  // We need to map the active element ids into a
  // contiguous range.
  {
//     active_elem_iterator       elem_it (mesh.elements_begin());
//     const active_elem_iterator elem_end(mesh.elements_end());

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 

    unsigned int el_num = 0;

    for (; elem_it != elem_end; ++elem_it)
      {
	assert ((*elem_it)->id() < forward_map.size());
	
	forward_map[(*elem_it)->id()] = el_num++;
      }

    assert (el_num == n_active_elem);
   }
  
  
  // build the graph in CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors
  {
    std::vector<const Elem*> neighbors_offspring;
    
//     const_active_elem_iterator       elem_it (mesh.const_elements_begin());
//     const const_active_elem_iterator elem_end(mesh.const_elements_end());

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 

    // This will be exact when there is no refinement and all the
    // elements are of the same type.
    adjncy.reserve (n_active_elem*(*elem_it)->n_neighbors());
    
    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;

	assert (elem->id() < forward_map.size());
	assert (forward_map[elem->id()] != libMesh::invalid_uint);
	
	// maybe there is a better weight?
	// The weight is used to define what a balanced graph is
	vwgt[forward_map[elem->id()]] = elem->n_nodes(); 

	// The beginning of the adjacency array for this elem
	xadj.push_back(adjncy.size());

	// Loop over the element's neighbors.  An element
	// adjacency corresponds to a face neighbor
	for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
	  {
	    const Elem* neighbor = elem->neighbor(ms);
	    
	    if (neighbor != NULL)
	      {  
		// If the neighbor is active treat it
		// as a connection
		if (neighbor->active())
		  {
		    assert (neighbor->id() < forward_map.size());
		    assert (forward_map[neighbor->id()] != libMesh::invalid_uint);

		    adjncy.push_back (forward_map[neighbor->id()]);
		  }
  
#ifdef ENABLE_AMR
		
		// Otherwise we need to find all of the
		// neighbor's children that are connected to
		// us and add them
		else
		  {
		    // The side of the neighbor to which
		    // we are connected
		    const unsigned int ns =
		      neighbor->which_neighbor_am_i (elem);
		    
		    // Get all the active children (& grandchildren, etc...)
		    // of the neighbor.
		    neighbor->active_family_tree (neighbors_offspring);
		    
		    // Get all the neighbor's children that
		    // live on that side and are thus connected
		    // to us
		    for (unsigned int nc=0; nc<neighbors_offspring.size(); nc++)
		      {
			const Elem* child =
			  neighbors_offspring[nc];
			
			// This does not assume a level-1 mesh.
			// Note that since children have sides numbered
			// coincident with the parent then this is a sufficient test.
			if (child->neighbor(ns) == elem)
			  {
			    assert (child->active());
			    assert (child->id() < forward_map.size());
			    assert (forward_map[child->id()] != libMesh::invalid_uint);
			
			    adjncy.push_back (forward_map[child->id()]);
			  }
		      }
		  }

#endif /* ifdef ENABLE_AMR */

	      }
	  }
      }
    
    // The end of the adjacency array for the last elem
    xadj.push_back(adjncy.size());
    
  } // done building the graph


  // Select which type of partitioning to create

  // Use recursive if the number of partitions is less than or equal to 8
  if (n_pieces <= 8)
    Metis::METIS_PartGraphRecursive(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
				    &wgtflag, &numflag, &nparts, &options[0],
				    &edgecut, &part[0]);

  // Otherwise  use kway
  else
    Metis::METIS_PartGraphKway(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
			       &wgtflag, &numflag, &nparts, &options[0],
			       &edgecut, &part[0]);
  
  
  // Assign the returned processor ids.  The part array contains
  // the processor id for each active element, but in terms of
  // the contiguous indexing we defined above
  {
 //    active_elem_iterator       elem_it (mesh.elements_begin());
//     const active_elem_iterator elem_end(mesh.elements_end());

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 

    for (; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;

	assert (elem->id() < forward_map.size());
	assert (forward_map[elem->id()] != libMesh::invalid_uint);
	
	elem->set_processor_id() =
	  static_cast<short int>(part[forward_map[elem->id()]]);	
      }
  }

  STOP_LOG("partition()", "MetisPartitioner");
  
#endif
  
}
