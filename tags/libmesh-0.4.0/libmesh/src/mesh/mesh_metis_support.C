// $Id: mesh_metis_support.C,v 1.15 2003-05-29 04:29:16 benkirk Exp $

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



// C++ includes
#include <math.h>

// Local includes
#include "mesh_base.h"


#ifdef HAVE_METIS

# include "libmesh.h"
# include "elem.h"
# include "mesh_logging.h"

  namespace Metis {
    extern "C" {
#     include "metis.h"
    }
  }

#endif


void MeshBase::metis_partition(const unsigned int n_sbdmns_in,
			       const std::string& type)
{
  const unsigned int n_active_elem = this->n_active_elem();
  const unsigned int n_sbdmns      = std::min (n_sbdmns_in, n_active_elem);
  
  this->set_n_subdomains() = n_sbdmns;

  // check for easy return
  if (n_sbdmns == 1)
    {
      elem_iterator       elem_it (this->elements_begin());
      const elem_iterator elem_end(this->elements_end());

      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() = 0;
      
      return;
    }

  
#ifndef HAVE_METIS

  
  
  std::cerr << "ERROR:  Metis not detected during configuration!" << std::endl
	    << "        Using space-filling curves instead."      << std::endl;

  here();
  
  this->sfc_partition(n_sbdmns);


  
#else

  

  assert (this->mesh_dimension() != 1);

  START_LOG("metis_partition()", "MeshBase");

  // build the graph
  // the forward_map maps the active element id
  // into a contiguous block of indices for Metis
  std::vector<unsigned int> forward_map (this->n_elem(),
					 static_cast<unsigned int>(-1));
  
  std::vector<int>          xadj;
  std::vector<int>          adjncy;
  std::vector<int>          options(5);
  std::vector<int>          vwgt(n_active_elem);
  std::vector<int>          part(n_active_elem);
  
  int
    n = static_cast<int>(n_active_elem),  // number of "nodes" (elements)
                                          //   in the graph
    wgtflag = 2,                          // weights on vertices only
    numflag = 0,                          // C-style 0-based numbering
    nparts  = static_cast<int>(n_sbdmns), // number of subdomains to create
    edgecut = 0;                          // the numbers of edges cut by the
                                          //   partition

  // Set the options
  options[0] = 0; // use default options


  // Metis will only consider the active elements.
  // We need to map the active element ids into a
  // contiguous range.
  {
    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());

    unsigned int el_num = 0;

    for (; elem_it != elem_end; ++elem_it)
      {
	assert ((*elem_it)->id() < forward_map.size());
	
	forward_map[(*elem_it)->id()] = el_num;
	el_num++;
      }

    assert (el_num == n_active_elem);
   }
  
  
  // build the graph in CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors
  {
    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());
    
    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;

	assert (elem->id() < forward_map.size());
	assert (forward_map[elem->id()] !=
		static_cast<unsigned int>(-1));
	
	// maybe there is a better weight?
	vwgt[forward_map[elem->id()]]
	     = elem->n_nodes(); 

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
		    assert (forward_map[neighbor->id()] !=
			    static_cast<unsigned int>(-1));

		    adjncy.push_back (forward_map[neighbor->id()]);

		  }
		
		// Otherwise we need to find all of the
		// neighbor's children that are connected to
		// us and add them
		else
		  {
		    // The side of the neighbor to which
		    // we are connected
		    const unsigned int ns =
		      neighbor->which_neighbor_am_i (elem);
		    
		    // Get all the neighbor's children that
		    // live on that side and are thus connected
		    // to us
		    for (unsigned int nc=0;
			 nc<neighbor->n_children_per_side(ns); nc++)
		      {
			const Elem* child =
			  neighbor->child (neighbor->side_children_matrix(ns,nc));
			
			// This assumes a level-1 mesh
			assert (child->active());
			assert (child->id() < forward_map.size());
			assert (forward_map[child->id()] !=
				static_cast<unsigned int>(-1));
			
			adjncy.push_back (forward_map[child->id()]);
		      }
		  }
	      }
	  }
      }
    
    // The end of the adjacency array for this elem
    xadj.push_back(adjncy.size());
    
  } // done building the graph


  // Select which type of partitioning to create

  // Use recursive if specified or if the number of partitions is
  // less than or equal to 8
  if ((type == "recursive") ||
      (n_sbdmns <= 8))
    Metis::METIS_PartGraphRecursive(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
				    &wgtflag, &numflag, &nparts, &options[0],
				    &edgecut, &part[0]);

  // Maybe use kway if specified
  else if (type == "kway")
    Metis::METIS_PartGraphKway(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
			       &wgtflag, &numflag, &nparts, &options[0],
			       &edgecut, &part[0]);
  

  // Otherwise throw an error message.
  else
    {
      std::cerr << " ERROR:  valid options:" << std::endl
		<< "   \"recursive\" " << std::endl
		<< "   \"kway\"  "     << std::endl
		<< " Using space-filling curves instead." << std::endl;

      STOP_LOG("metis_partition()", "MeshBase");
      sfc_partition(n_sbdmns);
      return;
    }

  
  // Assign the returned processor ids
  {
    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());

    for (; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;

	assert (elem->id() < forward_map.size());
	assert (forward_map[elem->id()] !=
		static_cast<unsigned int>(-1));
	
	elem->set_subdomain_id() =
	  elem->set_processor_id() =
	  static_cast<short int>(part[forward_map[elem->id()]]);
	
      }
  }

  STOP_LOG("metis_partition()", "MeshBase");

#endif
}
