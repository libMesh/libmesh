// $Id: partitioner.C,v 1.2 2003-08-26 22:58:45 jwpeterson Exp $

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
#include "partitioner.h"
#include "mesh_base.h"



// ------------------------------------------------------------
// Partitioner implementation
void Partitioner::single_partition ()
{
  // Loop over all the elements and assign them to processor 0.
  elem_iterator       elem_it (_mesh.elements_begin());
  const elem_iterator elem_end(_mesh.elements_end());
      
  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->set_processor_id() = 0;
}
