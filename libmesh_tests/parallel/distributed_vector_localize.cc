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
#include <distributed_vector.h>
#include <vector>



void unit_test ()
{
  unsigned int block_size  = 10;
  unsigned int local_size  = block_size + libMesh::processor_id();  // a different size on
  unsigned int global_size = 0;                                     // each processor.

  for (unsigned int p=0; p<libMesh::n_processors(); p++)
    global_size += (block_size + libMesh::processor_id());  

  {
    DistributedVector<Real> 
      v (global_size, local_size),
      l (global_size, global_size);
    
    const unsigned int 
      first = v.first_local_index(),
      last  = v.last_local_index();

    for (unsigned int n=first; n != last; n++)
      v.set (n, static_cast<Real>(n));
    v.close();
    here();
    v.localize(l);
    here();
    std::cout << " l.size() = " << l.size() 
	      << "v.size() = " << v.size()<< std::endl;
    std::cout << "v=[ ";
    for (unsigned int i=first; i<last; i++)
      std::cout << v(i) << " ";
    std::cout << "]" << std::endl;

    std::cout << "l=[ ";
    for (unsigned int i=0; i<l.size(); i++)
      std::cout << l(i) << " ";
    std::cout << "]" << std::endl;

  }
}
