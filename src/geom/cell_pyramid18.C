// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_pyramid18.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/fe_reference_element_traits.h"
#include "libmesh/face_tri7.h"
#include "libmesh/face_quad9.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid18 class static member initializations
const int Pyramid18::num_nodes;
const int Pyramid18::nodes_per_side;
const int Pyramid18::nodes_per_edge;

// ------------------------------------------------------------
// Pyramid18 class member functions

bool Pyramid18::is_vertex(const unsigned int i) const
{
  if (i < 5)
    return true;
  return false;
}



bool Pyramid18::is_edge(const unsigned int i) const
{
  if (i < 5)
    return false;
  if (i > 12)
    return false;
  return true;
}



bool Pyramid18::is_face(const unsigned int i) const
{
  if (i > 12)
    return true;
  return false;
}



bool Pyramid18::is_node_on_side(const unsigned int n,
                                const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  const auto count = side_node_count_or_zero(this->type(), s);
  for (unsigned int i = 0; i != count; ++i)
    if (this->local_side_node(s, i) == n)
      return true;
  return false;
}

std::vector<unsigned>
Pyramid18::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  const auto count = side_node_count_or_zero(this->type(), s);
  std::vector<unsigned> nodes(count);
  for (unsigned int i = 0; i != count; ++i)
    nodes[i] = this->local_side_node(s, i);
  return nodes;
}

std::vector<unsigned>
Pyramid18::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  const auto count = edge_node_count_or_zero(this->type(), e);
  std::vector<unsigned> nodes(count);
  for (unsigned int i = 0; i != count; ++i)
    nodes[i] = this->local_edge_node(e, i);
  return nodes;
}

bool Pyramid18::is_node_on_edge(const unsigned int n,
                                const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  const auto count = edge_node_count_or_zero(this->type(), e);
  for (unsigned int i = 0; i != count; ++i)
    if (this->local_edge_node(e, i) == n)
      return true;
  return false;
}



bool Pyramid18::has_affine_map() const
{
  // TODO: If the base is a parallelogram and all the triangular faces are planar,
  // the map should be linear, but I need to test this theory...
  return false;
}



Order Pyramid18::default_order() const
{
  return THIRD;
}



dof_id_type Pyramid18::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      return this->compute_key (this->node_id(s+14));

    case 4:  // the quad face at z=0
      return this->compute_key (this->node_id(13));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int Pyramid18::local_side_node(unsigned int side,
                                        unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 9 nodes per side.
  libmesh_assert_less(side_node, Pyramid18::nodes_per_side);

  // Some sides have 7 nodes.
  libmesh_assert(side == 4 || side_node < 7);

  unsigned int node = invalid_uint;
  libmesh_error_msg_if(!try_local_side_node(this->type(), side, side_node, node),
                       "Pyramid18::local_side_node(): unsupported shared side-node lookup");
  return node;
}



unsigned int Pyramid18::local_edge_node(unsigned int edge,
                                        unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, Pyramid18::nodes_per_edge);

  unsigned int node = invalid_uint;
  libmesh_error_msg_if(!try_local_edge_node(this->type(), edge, edge_node, node),
                       "Pyramid18::local_edge_node(): unsupported shared edge-node lookup");
  return node;
}



std::unique_ptr<Elem> Pyramid18::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        face = std::make_unique<Tri7>();
        break;
      }
    case 4: // the quad face at z=0
      {
        face = std::make_unique<Quad9>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(this->local_side_node(i, n)));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}



void Pyramid18::build_side_ptr (std::unique_ptr<Elem> & side,
                                const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        if (!side.get() || side->type() != TRI7)
          {
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }
    case 4:  // the quad face at z=0
      {
        if (!side.get() || side->type() != QUAD9)
          {
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->inherit_data_from(*this);

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n, this->node_ptr(this->local_side_node(i, n)));
}



std::unique_ptr<Elem> Pyramid18::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Pyramid18>(i);
}



void Pyramid18::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Pyramid18>(edge, i, EDGE3);
}



void Pyramid18::connectivity(const unsigned int libmesh_dbg_var(sc),
                             const IOPackage iop,
                             std::vector<dof_id_type> & /*conn*/) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        // TODO
        libmesh_not_implemented();
      }

    case VTK:
      {
        // TODO
        libmesh_not_implemented();
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



unsigned int Pyramid18::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
      return 2;

    case 13:
      return 4;

    case 14:
    case 15:
    case 16:
    case 17:
      return 3;

    default:
      libmesh_error_msg("Invalid node n = " << n);
    }
}


unsigned short int Pyramid18::second_order_adjacent_vertex (const unsigned int n,
                                                            const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
      {
        libmesh_assert_less (v, 2);

        // This is the analog of the static, const arrays
        // {Hex,Prism,Tet10}::_second_order_adjacent_vertices
        // defined in the respective source files... possibly treat
        // this similarly once the Pyramid13 has been added?
        constexpr unsigned short node_list[8][2] =
          {
            {0,1},
            {1,2},
            {2,3},
            {0,3},
            {0,4},
            {1,4},
            {2,4},
            {3,4}
          };

        return node_list[n-5][v];
      }

      // mid-face node on bottom
    case 13:
      {
        libmesh_assert_less (v, 4);

        // The vertex nodes surrounding node 13 are 0, 1, 2, and 3.
        // Thus, the v'th node is simply = v.
        return cast_int<unsigned short>(v);
      }

      // mid-face nodes on triangles
    case 14:
    case 15:
    case 16:
    case 17:
      {
        libmesh_assert_less (v, 3);

        constexpr unsigned short node_list[4][3] =
          {
            {0,1,4},
            {1,2,4},
            {2,3,4},
            {0,3,4}
          };

        return node_list[n-14][v];
      }

    default:
      libmesh_error_msg("Invalid n = " << n);

    }
}



void Pyramid18::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(5,6,7,8);
      swap4nodes(9,10,11,12);
      swap4nodes(14,15,16,17);
      swap4neighbors(0,1,2,3);
    }
}


void Pyramid18::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2nodes(6,8);
  swap2nodes(9,10);
  swap2nodes(11,12);
  swap2nodes(15,17);
  swap2neighbors(1,3);
  swap2boundarysides(1,3,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
  swap2boundaryedges(4,5,boundary_info);
  swap2boundaryedges(6,7,boundary_info);
}


unsigned int Pyramid18::center_node_on_side(const unsigned short side) const
{
  libmesh_assert_less (side, Pyramid18::num_sides);
  return side == 4 ? 13 : side+14;
}


ElemType Pyramid18::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s < 4)
    return TRI7;
  return QUAD9;
}


} // namespace libMesh
