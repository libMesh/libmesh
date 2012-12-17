// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_NODE_RANGE_H
#define LIBMESH_NODE_RANGE_H

// Local includes
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/stored_range.h"

// C++ includes

namespace libMesh
{

typedef StoredRange<MeshBase::node_iterator,             Node*>      NodeRange;
typedef StoredRange<MeshBase::const_node_iterator, const Node*> ConstNodeRange;

} // namespace libMesh

#endif // LIBMESH_NODE_RANGE_H
