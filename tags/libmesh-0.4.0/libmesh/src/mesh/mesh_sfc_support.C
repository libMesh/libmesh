// $Id: mesh_sfc_support.C,v 1.1 2003-05-29 04:29:16 benkirk Exp $

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
#include "elem.h"
#include "mesh_logging.h"


#ifdef HAVE_SFCURVES

  namespace Sfc {
    extern "C" {
#     include "sfcurves.h"
    }
  }
#endif




// ------------------------------------------------------------
// MeshBase class member functions
void MeshBase::sfc_partition(const unsigned int n_sbdmns_in,
			     const std::string& type)
{

  const unsigned int n_active_elem = this->n_active_elem();
  const unsigned int n_sbdmns      = std::min (n_sbdmns_in, n_active_elem);
  
  // Set the number of subdomains
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

  
  
  // won't work without SFC library!
#ifndef HAVE_SFCURVES
  
  std::cerr << "ERROR: Not compiled with space-filling curve"      << std::endl
	    << "       support. Using linear partitioning instead" << std::endl
	    << "       This partitioning could be arbitrarily bad!"<< std::endl
	    << std::endl;

  here();

  // Create a simple linear partitioning
  {
    const unsigned int blksize = n_active_elem/n_sbdmns; 
    unsigned int e = 0;

    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());

    for ( ; elem_it != elem_end; ++elem_it)
      {
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() = 
	  e/blksize;
	
	e++;
      }
  }

  
#else

  
  START_LOG("sfc_partition()", "MeshBase");

  // the forward_map maps the active element id
  // into a contiguous block of indices for Metis
  std::vector<unsigned int> forward_map (this->n_elem(),
					 static_cast<unsigned int>(-1));

  // the reverse_map maps the contiguous ids back
  // to active elements
  std::vector<Elem*> reverse_map (n_active_elem, NULL);
  
  int size = static_cast<int>(n_active_elem);
  std::vector<double> x     (size);
  std::vector<double> y     (size);
  std::vector<double> z     (size);
  std::vector<int>    table (size);


  // We need to map the active element ids into a
  // contiguous range.
  {
    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());

    unsigned int el_num = 0;

    for (; elem_it != elem_end; ++elem_it)
      {
	assert ((*elem_it)->id() < forward_map.size());
	assert (el_num           < reverse_map.size());
	
	forward_map[(*elem_it)->id()] = el_num;
	reverse_map[el_num]           = *elem_it;
	el_num++;
      }
    assert (el_num == n_active_elem);
   }


  // Get the centroid for each active element
  {
    active_elem_iterator       elem_it (this->elements_begin());
    const active_elem_iterator elem_end(this->elements_end());

    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;
	
	assert (elem->id() < forward_map.size());

	const Point p = elem->centroid();

	x[forward_map[elem->id()]] = p(0);
	y[forward_map[elem->id()]] = p(1);
	z[forward_map[elem->id()]] = p(2);
      }
  }


  // build the space-filling curve
  if (type == "hilbert")
    Sfc::hilbert(&x[0], &y[0], &z[0], &size, &table[0]);
  
  else if (type == "morton")
    Sfc::morton(&x[0], &y[0], &z[0], &size, &table[0]);
  
  else
    error();

  
  // Assign the partitioning to the active elements
  {
    const unsigned int blksize = n_active_elem/n_sbdmns; 

    for (unsigned int i=0; i<n_active_elem; i++)
      {
	assert (static_cast<unsigned int>(table[i]-1) < reverse_map.size());
		
	Elem* elem = reverse_map[table[i]-1];

	elem->set_subdomain_id() = 
	  elem->set_processor_id() =
	  i/blksize;
      }
  }
  
  STOP_LOG("sfc_partition()", "MeshBase");
  
#endif
}




