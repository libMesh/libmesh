
// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_PARALLEL_ELEM_H
#define LIBMESH_PARALLEL_ELEM_H

// Local Includes
#include "libmesh/libmesh_config.h"

#include "libmesh/elem.h"
#include "libmesh/id_types.h"

namespace libMesh {
namespace Parallel {
// BufferType<> specializations to return a buffer datatype
// to handle communication of Elems
template <>
struct BufferType<const Elem*> {
  typedef largest_id_type type;
};

template <>
struct BufferType<Elem> {
  typedef largest_id_type type;
};

} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_ELEM_H
