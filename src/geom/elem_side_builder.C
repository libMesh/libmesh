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



#include "libmesh/elem_side_builder.h"

namespace libMesh
{

Elem &
ElemSideBuilder::operator()(Elem & elem, const unsigned int s)
{
  libmesh_assert_less(s, elem.n_sides());
  const std::size_t type_index = static_cast<std::size_t>(elem.side_type(s));
  if (type_index >= _cached_elems.size())
    _cached_elems.resize(type_index + 1);
  std::unique_ptr<Elem> & side_elem = _cached_elems[type_index];
  elem.build_side_ptr(side_elem, s);
  return *side_elem;
}

const Elem &
ElemSideBuilder::operator()(const Elem & elem, const unsigned int s)
{
  return (*this)(const_cast<Elem &>(elem), s);
}

} // namespace libMesh
