// $Id: linear_partitioner.C,v 1.4 2003-09-02 18:02:44 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh.h"
#include "linear_partitioner.h"
#include "mesh_logging.h"



// ------------------------------------------------------------
// LinearPartitioner implementation
void LinearPartitioner::partition (const unsigned int n)
{
  assert (n > 0);

  // Check for an easy return
  if (n == 1)
    {
      this->single_partition ();
      return;
    }
  
  // Create a simple linear partitioning
  {
    START_LOG ("partition()", "LinearPartitioner");
    
    const unsigned int n_active_elem = _mesh.n_active_elem();
    const unsigned int blksize       = n_active_elem/n;
    
    unsigned int e = 0;
        
    active_elem_iterator       elem_it (_mesh.elements_begin());
    const active_elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      {
	(*elem_it)->set_processor_id() = e/blksize;
	
	e++;
      }
    
    STOP_LOG ("partition()", "LinearPartitioner");
  }
}
