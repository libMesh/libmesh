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

#ifndef LIBMESH_CELL_HEX8_IMPL_H
#define LIBMESH_CELL_HEX8_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hex8 class member functions

template <typename RealType>
bool Hex8Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Hex8Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Hex8Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Hex8Templ<RealType>::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Hex8Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

template <typename RealType>
bool Hex8Templ<RealType>::is_node_on_edge(const unsigned int n,
                           const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



template <typename RealType>
bool Hex8Templ<RealType>::has_affine_map() const
{
  // Make sure x-edge endpoints are affine
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(4)) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(7)))
    return false;
  // Make sure xz-faces are identical parallelograms
  v = this->point(4) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(3)))
    return false;
  // If all the above checks out, the map is affine
  return true;
}



template <typename RealType>
Order Hex8Templ<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Hex8Templ<RealType>::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Quad4,Hex8>>(this,i);

  else
    {
      std::unique_ptr<Elem> face = libmesh_make_unique<Quad4>();
      face->subdomain_id() = this->subdomain_id();

      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Hex8::side_nodes_map[i][n]);

      return face;
    }
}



template <typename RealType>
void Hex8Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->template simple_build_side_ptr<Hex8>(side, i, QUAD4);
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Hex8Templ<RealType>::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2,Hex8>>(this,i);
}



template <typename RealType>
void Hex8Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sc),
                                       const IOPackage iop,
                                       std::vector<dof_id_type> & conn) const
{
  libmesh_assert(this->_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(8);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        conn[4] = this->node_id(4)+1;
        conn[5] = this->node_id(5)+1;
        conn[6] = this->node_id(6)+1;
        conn[7] = this->node_id(7)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        conn[3] = this->node_id(3);
        conn[4] = this->node_id(4);
        conn[5] = this->node_id(5);
        conn[6] = this->node_id(6);
        conn[7] = this->node_id(7);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}

template <typename RealType>
RealType Hex8Templ<RealType>::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0), x1 = this->point(1), x2 = this->point(2), x3 = this->point(3),
    x4 = this->point(4), x5 = this->point(5), x6 = this->point(6), x7 = this->point(7);

  // Construct constant data vectors.  The notation is:
  // \vec{x}_{\xi}   = \vec{a1}*eta*zeta + \vec{b1}*eta + \vec{c1}*zeta + \vec{d1}
  // \vec{x}_{\eta}  = \vec{a2}*xi*zeta  + \vec{b2}*xi  + \vec{c2}*zeta + \vec{d2}
  // \vec{x}_{\zeta} = \vec{a3}*xi*eta   + \vec{b3}*xi  + \vec{c3}*eta  + \vec{d3}
  // but it turns out that a1, a2, and a3 are not needed for the volume calculation.

  // Build up the 6 unique vectors which make up dx/dxi, dx/deta, and dx/dzeta.
  Point q[6] =
    {
      /*b1*/  x0 - x1 + x2 - x3 + x4 - x5 + x6 - x7, /*=b2*/
      /*c1*/  x0 - x1 - x2 + x3 - x4 + x5 + x6 - x7, /*=b3*/
      /*d1*/ -x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7,
      /*c2*/  x0 + x1 - x2 - x3 - x4 - x5 + x6 + x7, /*=c3*/
      /*d2*/ -x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7,
      /*d3*/ -x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7
    };

  // We could check for a linear element, but it's probably faster to
  // just compute the result...
  return
    (triple_product(q[0], q[4], q[3]) +
     triple_product(q[2], q[0], q[1]) +
     triple_product(q[1], q[3], q[5])) / 192. +
    triple_product(q[2], q[4], q[5]) / 64.;
}

template <typename RealType>
BoundingBoxTempl<RealType>
Hex8Templ<RealType>::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}

} // namespace libMesh

#endif // LIBMESH_CELL_HEX8_IMPL_H
