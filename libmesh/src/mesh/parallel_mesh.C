// $Id: parallel_mesh.C,v 1.4 2007-09-21 19:55:55 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "elem.h"
#include "parallel_mesh.h"

// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (unsigned int d) :
  Mesh (d)
{
}


ParallelMesh::~ParallelMesh ()
{
}


void ParallelMesh::delete_nonlocal_elements()
{
  std::vector<bool> local_nodes(_nodes.size(), false);
  std::vector<bool> semilocal_elems(_elements.size(), false);

  const_element_iterator l_elem_it = this->local_elements_begin(),
                         l_end     = this->local_elements_end();
  for (; l_elem_it != l_end; ++l_elem_it)
    {
      const Elem *elem = *l_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
    }

  element_iterator nl_elem_it = this->not_local_elements_begin(),
                   nl_end     = this->not_local_elements_end();
  for (; nl_elem_it != nl_end; ++nl_elem_it)
    {
      Elem *elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            semilocal_elems[elem->id()] = true;
            break;
          }
      if (!semilocal_elems[elem->id()])
        {
          // delete_elem doesn't currently invalidate element
          // iterators... that had better not change
          this->delete_elem(elem);
        }
    }
//  We eventually want to compact the _nodes and _elems vectors
//  But not do a repartition or anything
//  do not this->prepare_for_use();
}

void ParallelMesh::restore_nonlocal_elements()
{
  // Someday...
  error();
}
