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
#include "libmesh/edge_edge4.h"
#include "libmesh/face_quad16.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// Quad16 class static member initialization
const int Quad16::num_nodes;
const int Quad16::nodes_per_side;


// Note: we can omit initialization of the third entry of each row because
// static variables are automatically zero-initialized.
// After the vertices the node locations start moving from Quad8/9

static constexpr Real _third = Real(1)/3;
const Real Quad16::_master_points[16][3] =
{
  //corners
  {-1, -1},
  {1, -1},
  {1, 1},
  {-1, 1},
  //internal edges
  {- _third, -1},
  {_third, -1},
  //right
  {1, - _third},
  {1, _third},
  //top
  {_third, 1},
  {- _third, 1},
  // left
  {-1, _third},
  {-1, -_third},
  //internal
  {- _third, - _third},
  {_third, - _third},
  {_third, _third},
  {- _third, _third}
};

const unsigned int Quad16::side_nodes_map[Quad16::num_sides][Quad16::nodes_per_side] =
  {
    {0,  1,  4,  5}, // Side 0
    {1,  2,  6,  7}, // Side 1
    {2,  3,  8,  9}, // Side 2
    {3,  0, 10, 11}  // Side 3
  };

// ------------------------------------------------------------
// Quad16 class member functions

bool Quad16::is_vertex(const unsigned int n) const
{
  if (n < 4)
    return true;
  return false;
}

bool Quad16::is_edge(const unsigned int n) const
{
  libmesh_assert_not_equal_to (n, invalid_uint);
  if (n < 4)
    return false;
  if (n > 11)
    return false;
  return true;
}

bool Quad16::is_face(const unsigned int n) const
{
  libmesh_assert_not_equal_to (n, invalid_uint);
  if (n > 11)
    return true;
  return false;
}

bool Quad16::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

unsigned int Quad16::local_side_node(unsigned int side,
                                   unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Quad16::nodes_per_side);

  return Quad16::side_nodes_map[side][side_node];
}

std::vector<unsigned>
Quad16::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Quad16::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
}

bool Quad16::has_affine_map() const
{
  Point v = this->point(3) - this->point(0);
  return (v.relative_fuzzy_equals(this->point(2) - this->point(1), affine_tol));
}


std::unique_ptr<Elem> Quad16::build_side_ptr (const unsigned int i)
{
  return this->simple_build_side_ptr<Edge4, Quad16>(i);
}



void Quad16::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->simple_build_side_ptr<Quad16>(side, i, EDGE4);
}



void Quad16::connectivity(const unsigned int libmesh_dbg_var(sf),
                         const IOPackage iop,
                         std::vector<dof_id_type> &) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);
  libmesh_error_msg("Unsupported IO package for Quad16 " << iop);
}


void Quad16::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(4,6,8,10);
      swap4nodes(5,7,9,11);
      swap4nodes(12,13,14,15);
      swap4neighbors(0,1,2,3);
    }
}


void Quad16::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2nodes(4,5);
  swap2nodes(6,11);
  swap2nodes(7,10);
  swap2nodes(8,9);
  swap2nodes(12,13);
  swap2nodes(14,15);
  swap2neighbors(1,3);
  swap2boundarysides(1,3,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
}

ElemType Quad16::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 4);
  return EDGE4;
}


} // namespace libMesh
