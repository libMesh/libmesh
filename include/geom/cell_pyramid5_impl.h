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

#ifndef LIBMESH_CELL_PYRAMID5_IMPL_H
#define LIBMESH_CELL_PYRAMID5_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Pyramid5 class member functions

template <typename RealType>
bool Pyramid5Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Pyramid5Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Pyramid5Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Pyramid5Templ<RealType>::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Pyramid5Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  auto trim = (s == 4) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

template <typename RealType>
bool Pyramid5Templ<RealType>::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



template <typename RealType>
bool Pyramid5Templ<RealType>::has_affine_map() const
{
  //  Point v = this->point(3) - this->point(0);
  //  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
  return false;
}



template <typename RealType>
Order Pyramid5Templ<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Pyramid5Templ<RealType>::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:
        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<Tri3Templ<Real>,Pyramid5>>(this,i);

        case 4:
          return libmesh_make_unique<Side<Quad4Templ<Real>,Pyramid5>>(this,i);

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Return value
      std::unique_ptr<Elem> face;

      switch (i)
        {
        case 0: // triangular face 1
        case 1: // triangular face 2
        case 2: // triangular face 3
        case 3: // triangular face 4
          {
            face = libmesh_make_unique<Tri3Templ<Real>>();
            break;
          }
        case 4: // the quad face at z=0
          {
            face = libmesh_make_unique<Quad4Templ<Real>>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Pyramid5Templ<RealType>::side_nodes_map[i][n]);

      return face;
    }
}



template <typename RealType>
void Pyramid5Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  this->side_ptr(side, i);
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Pyramid5Templ<RealType>::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2Templ<Real>,Pyramid5>>(this,i);
}



template <typename RealType>
void Pyramid5Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sc),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert(this->_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        conn[4] = this->node_id(4)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(4)+1;
        conn[7] = this->node_id(4)+1;
        return;
      }

    case VTK:
      {
        conn.resize(5);
        conn[0] = this->node_id(3);
        conn[1] = this->node_id(2);
        conn[2] = this->node_id(1);
        conn[3] = this->node_id(0);
        conn[4] = this->node_id(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


template <typename RealType>
RealType Pyramid5Templ<RealType>::volume () const
{
  // The pyramid with a bilinear base has volume given by the
  // formula in: "Calculation of the Volume of a General Hexahedron
  // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
  Point
    x0 = this->point(0), x1 = this->point(1), x2 = this->point(2),
    x3 = this->point(3), x4 = this->point(4);

  // Construct various edge and diagonal vectors.
  Point v40 = x0 - x4;
  Point v13 = x3 - x1;
  Point v02 = x2 - x0;
  Point v03 = x3 - x0;
  Point v01 = x1 - x0;

  // Finally, ready to return the volume!
  return
    triple_product(v40, v13, v02) / 6. +
    triple_product(v02, v01, v03) / 12.;
}

template <typename RealType>
BoundingBoxTempl<RealType>
Pyramid5Templ<RealType>::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}

} // namespace libMesh

#endif // LIBMESH_CELL_PYRAMID5_IMPL_H
