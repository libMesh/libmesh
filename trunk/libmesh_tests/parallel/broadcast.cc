// $Id: allgather.cc 2236 2007-10-24 18:27:40Z benkirk $
 
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

#include <parallel.h>
#include <vector>

void unit_test ()
{
  std::vector<unsigned int> src(3), dest(3);

  src[0]=0;
  src[1]=1;
  src[2]=2;

  if (libMesh::processor_id() == 0)
    dest = src;

  Parallel::broadcast(dest);

  if (src != dest)
    {
      std::cerr << "Test failed!" << std::endl;
      error();
    }

  std::cout << "Test succeeded!" << std::endl;
}
