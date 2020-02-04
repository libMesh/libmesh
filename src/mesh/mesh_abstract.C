// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute and/or
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

#include "libmesh/mesh_abstract.h"

#include <map>
#include <set>

namespace libMesh
{

MeshAbstract::MeshAbstract(const Parallel::Communicator & comm_in) :
    ParallelObject (comm_in) {}

unsigned int MeshAbstract::mesh_dimension() const
{
  if (!_elem_dims.empty())
    return cast_int<unsigned int>(*_elem_dims.rbegin());
  return 0;
}

}
