// $Id: vector_union.cc 2236 2007-10-24 18:27:40Z benkirk $
 
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

#include <libmesh.h>
#include <parallel.h>
#include <vector>



void unit_test ()
{
  const unsigned int root_id = 0;

  std::vector<unsigned int> v(3);

  for (unsigned int i=0; i<v.size(); i++)
    v[i] = v.size()*libMesh::processor_id() + i;

  if (libMesh::processor_id() == 0)
    v.push_back(libMesh::invalid_uint);

  // v should only expand on the root processor!
  const unsigned int old_size = v.size();

  Parallel::gather (root_id, v);

  if (libMesh::processor_id() != root_id)
    assert (v.size() == old_size);

  else // Global vector only on processor root_id!
    {
      std::cout << "v=[ ";
      for (unsigned int i=0; i<v.size(); i++)
	std::cout << v[i] << " ";

      std::cout << "]" << std::endl;
    }
}
