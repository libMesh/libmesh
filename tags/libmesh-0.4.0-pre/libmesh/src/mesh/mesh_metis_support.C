// $Id: mesh_metis_support.C,v 1.11 2003-03-04 15:31:24 benkirk Exp $

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
#include "libmesh.h"
#include "elem.h"


#ifdef HAVE_METIS
namespace Metis {
extern "C" {
# include "metis.h"
}
}
#endif


void MeshBase::metis_partition(const unsigned int n_sbdmns,
			       const std::string& type)
{
#ifndef HAVE_METIS
  
  std::cerr << "ERROR:  Metis not detected during configuration!" << std::endl
	    << "        Using space-filling curves instead." << std::endl;

  sfc_partition(n_sbdmns);
  
  return;
  
#else

  assert (n_sbdmns <= n_elem());
  
  set_n_subdomains() = n_sbdmns;
  set_n_processors() = n_sbdmns;

  // check for easy return
  if (n_sbdmns == 1)
    {
      for (unsigned int e=0; e<n_elem(); e++)
	elem(e)->set_subdomain_id() = 
	  elem(e)->set_processor_id() = 0;
      
      return;
    }
  
  
  assert (_dim != 1);

  START_LOG("metis_partition()", "MeshBase");

  // new way, build the graph
  std::vector<int> xadj;
  std::vector<int> adjncy;
  std::vector<int> options(5);
  std::vector<int> vwgt(n_elem());
  std::vector<int> part(n_elem());  
  
  int
    n = static_cast<int>(n_elem()),      // number of "nodes" (elements)
                                         //  in the graph
    wgtflag = 2,                         // weights on vertices only
    numflag = 0,                         // C-style 0-based numbering
    nparts = static_cast<int>(n_sbdmns), // number of subdomains to create
    edgecut = 0;                         // the numbers of edges cut by the
                                         //  partition

  options[0] = 0; // use default options

  // build the graph in CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors
  {
    std::map<const Elem*, int> elem_numbers;
    
    for (unsigned int e=0; e<n_elem(); e++)
      {
        vwgt[e] = elem(e)->n_nodes(); // maybe there is a better weight? 
        elem_numbers[elem(e)] = static_cast<int>(e);
      }

    bool found_a_neighbor = false;
    
    for (unsigned int e=0; e<n_elem(); e++)
      {
        xadj.push_back(adjncy.size());
        for (unsigned int s=0; s<elem(e)->n_sides(); s++)
          {
            const Elem* neighbor = elem(e)->neighbor(s);
            if (neighbor != NULL)
              {
                found_a_neighbor = true;
                adjncy.push_back(elem_numbers[neighbor]);
              }
          }
      }
    xadj.push_back(adjncy.size());

    // If we didn't find _any_ neighbors it is a safe bet the
    // user didn't call find_neighbors().  Print an informative
    // message and use a space-filling curve instead.
    if (!found_a_neighbor)
      {
	std::cerr << "ERROR: something is amiss..." << std::endl
		  << " I couldn't find any elements with neighbors." << std::endl
		  << " Did you forget to call find_neighbors()" << std::endl
		  << " BEFORE calling the graph partitioner?" << std::endl << std::endl
		  << " I'll use a space-filling curves instead." << std::endl;
	
	STOP_LOG("metis_partition()", "MeshBase");
	sfc_partition(n_sbdmns);
	return;
      }
  } // done with the map.


  if ((type == "recursive") ||
      (n_sbdmns <= 8))
    Metis::METIS_PartGraphRecursive(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
				    &wgtflag, &numflag, &nparts, &options[0],
				    &edgecut, &part[0]);

  else if (type == "kway")
    Metis::METIS_PartGraphKway(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
			       &wgtflag, &numflag, &nparts, &options[0],
			       &edgecut, &part[0]);
  

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


  for (unsigned int e=0; e<n_elem(); e++)
    elem(e)->set_subdomain_id() = 
      elem(e)->set_processor_id() = 
      static_cast<short int>(part[e]);

  STOP_LOG("metis_partition()", "MeshBase");

#endif
}
