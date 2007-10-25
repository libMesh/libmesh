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

#include <libmesh.h>
#include <parallel.h>
#include <vector>



void unit_test ()
{
  std::vector<unsigned int> v(3);

  for (unsigned int i=0; i<v.size(); i++)
    v[i] = v.size()*libMesh::processor_id() + i;

  if (libMesh::processor_id() == 0)
    v.push_back(libMesh::invalid_uint);

  Parallel::allgather (v);

  std::cout << "v=[ ";
  for (unsigned int i=0; i<v.size(); i++)
    std::cout << v[i] << " ";

  std::cout << "]" << std::endl;
}
