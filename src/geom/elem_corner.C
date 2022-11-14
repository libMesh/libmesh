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
#include "libmesh/elem_corner.h"

#if LIBMESH_DIM > 1

namespace libMesh
{

#if LIBMESH_DIM > 2
bool ElemCorner::at_edge(const Elem & elem, const unsigned short e) const
{
  libmesh_assert_less(e, elem.n_edges());
  if (!at_edge())
    return false;
  const unsigned int * nodes_on_edge_ptr = elem.nodes_on_edge_ptr(e);
  return at_edge(nodes_on_edge_ptr[0], nodes_on_edge_ptr[1]);
}

std::unique_ptr<const Elem> ElemCorner::build_edge(const Elem & elem) const
{
  libmesh_assert(at_edge());
  libmesh_assert_less(first, elem.n_vertices());
  libmesh_assert_less(second, elem.n_vertices());

  for (const auto e : elem.edge_index_range())
    if (elem.is_node_on_edge(first, e) && elem.is_node_on_edge(second, e))
      return elem.build_edge_ptr(e);

  libmesh_error_msg("Element does not contain vertices in ElemCorner");
}
#endif

std::string ElemCorner::print() const
{
  std::stringstream oss;
  if (at_vertex())
    oss << "at vertex " << vertex();
#if LIBMESH_DIM > 2
  else if (at_edge())
    oss << "at edge with vertices " << first << " and " << second;
#endif
  else
    oss << "not at corner";
  return oss.str();
}

bool ElemCorner::is_valid(const Elem & elem, const Point & point, const Real tol) const
{
  libmesh_assert(first == Elem::invalid_vertex || first < elem.n_vertices());
  libmesh_assert(second == Elem::invalid_vertex || second < elem.n_vertices());

  if (at_vertex())
    return vertex_point(elem).absolute_fuzzy_equals(point, tol);
#if LIBMESH_DIM > 2
  if (elem.dim() == 3 && at_edge())
    return build_edge(elem)->contains_point(point, tol);
#endif

  return true;
}

} // namespace libMesh

std::ostream &
operator<<(std::ostream & os, const libMesh::ElemCorner & elem_corner)
{
  os << elem_corner.print();
  return os;
}

#endif // LIBMESH_DIM > 1
