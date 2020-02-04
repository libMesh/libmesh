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



#ifndef LIBMESH_MESH_INSERTER_ITERATOR_H
#define LIBMESH_MESH_INSERTER_ITERATOR_H

// Local includes
#include "libmesh/mesh_base.h"

// C++ includes
#include <iterator>

namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
template <typename> class NodeTempl;
typedef NodeTempl<Real> Node;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;

/**
 * A class for templated methods that expect output iterator
 * arguments, which adds objects to the Mesh.
 * Although any mesh_inserter_iterator can add any object, we
 * template it around object type so that type inference and
 * iterator_traits will work.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief An output iterator for use with packed_range functions.
 */
template <typename T>
struct mesh_inserter_iterator
  : std::iterator<std::output_iterator_tag, T>
{
  mesh_inserter_iterator (MeshBase & m) : mesh(m) {}

  void operator=(Elem * e) { mesh.add_elem(e); }

  void operator=(Node * n) { mesh.insert_node(n); }

  void operator=(Point * p) { mesh.add_point(*p); }

  mesh_inserter_iterator & operator++() {
    return *this;
  }

  mesh_inserter_iterator operator++(int) {
    return mesh_inserter_iterator(*this);
  }

  // We don't return a reference-to-T here because we don't want to
  // construct one or have any of its methods called.  We just want
  // to allow the returned object to be able to do mesh insertions
  // with operator=().
  mesh_inserter_iterator & operator*() { return *this; }
private:

  MeshBase & mesh;
};

} // namespace libMesh


#endif // LIBMESH_MESH_INSERTER_ITERATOR_H
