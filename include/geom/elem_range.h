// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ELEM_RANGE_H
#define LIBMESH_ELEM_RANGE_H

// Local includes
#include "libmesh/mesh_base.h"
#include "libmesh/stored_range.h"

namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;

typedef StoredRange<MeshBase::element_iterator,             Elem *>      ElemRange;
typedef StoredRange<MeshBase::const_element_iterator, const Elem *> ConstElemRange;

template <typename RealType>
using ElemRangeTempl = StoredRange<typename MeshBaseTempl<RealType>::element_iterator, ElemTempl<RealType> *>;
template <typename RealType>
using ConstElemRangeTempl = StoredRange<typename MeshBaseTempl<RealType>::const_element_iterator, const ElemTempl<RealType> *>;

} // namespace libMesh

#endif // LIBMESH_ELEM_RANGE_H
