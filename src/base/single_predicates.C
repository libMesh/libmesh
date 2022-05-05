// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/single_predicates.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/distributed_mesh.h"

namespace libMesh
{

namespace Predicates
{

// The bid predicate returns true if has_boundary_id(node, id) returns true.
template <typename T>
bool
bid<T>::operator()(const T & it) const
{
  return _bndry_info.has_boundary_id(*it, _bid);
}


// The bnd predicate returns true if n_boundary_ids(node) > 0.
template <typename T>
bool
bnd<T>::operator()(const T & it) const
{
  return (_bndry_info.n_boundary_ids(*it) > 0);
}

template <typename T>
bool
evaluable<T>::operator()(const T & it) const
{
  return _dof_map.is_evaluable(**it, _var_num);
}

template <typename T>
bool
multi_evaluable<T>::operator()(const T & it) const
{
  for (const auto * const dof_map : _dof_maps)
  {
    libmesh_assert(dof_map);
    if (!dof_map->is_evaluable(**it))
      return false;
  }

  return true;
}


// Instantiate with the useful values of T
#define INSTANTIATE_NODAL_PREDICATES(IterType)                                               \
  template LIBMESH_EXPORT bool bid<IterType>::operator()(const IterType &) const;            \
  template LIBMESH_EXPORT bool bnd<IterType>::operator()(const IterType &) const;            \
  template LIBMESH_EXPORT bool evaluable<IterType>::operator()(const IterType &) const;      \
  template LIBMESH_EXPORT bool multi_evaluable<IterType>::operator()(const IterType &) const

#define INSTANTIATE_ELEM_PREDICATES(IterType)                                                \
  template LIBMESH_EXPORT bool evaluable<IterType>::operator()(const IterType &) const;      \
  template LIBMESH_EXPORT bool multi_evaluable<IterType>::operator()(const IterType &) const

// Handle commas in macro arguments
#define LIBMESH_COMMA ,

// This should probably be replaced with
// FooMesh::element_iterator_imp, etc
INSTANTIATE_ELEM_PREDICATES(std::vector<Elem *>::iterator);
INSTANTIATE_ELEM_PREDICATES(std::vector<Elem *>::const_iterator);
INSTANTIATE_NODAL_PREDICATES(std::vector<Node *>::iterator);
INSTANTIATE_NODAL_PREDICATES(std::vector<Node *>::const_iterator);
INSTANTIATE_ELEM_PREDICATES(DistributedMesh::dofobject_container<Elem>::veclike_iterator);
INSTANTIATE_ELEM_PREDICATES(DistributedMesh::dofobject_container<Elem>::const_veclike_iterator);
INSTANTIATE_NODAL_PREDICATES(DistributedMesh::dofobject_container<Node>::veclike_iterator);
INSTANTIATE_NODAL_PREDICATES(DistributedMesh::dofobject_container<Node>::const_veclike_iterator);


} // namespace Predicates

} // namespace libMesh
