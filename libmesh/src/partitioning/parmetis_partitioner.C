// $Id: parmetis_partitioner.C,v 1.1 2003-06-24 05:33:51 benkirk Exp $

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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "mesh_config.h"
#include "mesh.h"
#include "parmetis_partitioner.h"
#include "mesh_logging.h"

#ifdef HAVE_PARMETIS
  namespace Parmetis {
    extern "C" {
#     include "parmetis.h"
    }
  }
#else
#  include "metis_partitioner.h"
#endif


// ------------------------------------------------------------
// ParmetisPartitioner implementation
void ParmetisPartitioner::partition (const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      elem_iterator       elem_it (_mesh.elements_begin());
      const elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() =
	  0;
      
      return;
    }

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp(_mesh);

  mp.partition (n_sbdmns);
  
// What to do if the Metis library IS present
#else

  error();
  
  START_LOG("partition()", "ParmetisPartitioner");

  STOP_LOG("partition()", "ParmetisPartitioner");
  
#endif
  
}



void ParmetisPartitioner::repartition (const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      elem_iterator       elem_it (_mesh.elements_begin());
      const elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() =
	  0;
      
      return;
    }

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp(_mesh);

  mp.partition (n_sbdmns);
  
// What to do if the Metis library IS present
#else

  error();
  
  START_LOG("repartition()", "ParmetisPartitioner");

  STOP_LOG("repartition()", "ParmetisPartitioner");
  
#endif
  
}
