// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/mapvector.h"

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


// Instantiate with the useful values of T
#define INSTANTIATE_NODAL_PREDICATES(IterType)                          \
  template bool bid<IterType>::operator()(const IterType &) const;      \
  template bool bnd<IterType>::operator()(const IterType &) const;      \
  template bool evaluable<IterType>::operator()(const IterType &) const

#define INSTANTIATE_ELEM_PREDICATES(IterType)                           \
  template bool evaluable<IterType>::operator()(const IterType &) const

// Handle commas in macro arguments
#define LIBMESH_COMMA ,

// This should probably be replaced with
// FooMesh::element_iterator_imp, etc
INSTANTIATE_ELEM_PREDICATES(std::vector<Elem *>::iterator);
INSTANTIATE_ELEM_PREDICATES(std::vector<Elem *>::const_iterator);
INSTANTIATE_NODAL_PREDICATES(std::vector<Node *>::iterator);
INSTANTIATE_NODAL_PREDICATES(std::vector<Node *>::const_iterator);
INSTANTIATE_ELEM_PREDICATES(mapvector<Elem * LIBMESH_COMMA dof_id_type>::veclike_iterator);
INSTANTIATE_ELEM_PREDICATES(mapvector<Elem * LIBMESH_COMMA dof_id_type>::const_veclike_iterator);
INSTANTIATE_NODAL_PREDICATES(mapvector<Node * LIBMESH_COMMA dof_id_type>::veclike_iterator);
INSTANTIATE_NODAL_PREDICATES(mapvector<Node * LIBMESH_COMMA dof_id_type>::const_veclike_iterator);


} // namespace Predicates

} // namespace libMesh
