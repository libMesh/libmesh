// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/side.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri7.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace {
  constexpr libMesh::Real r18 = 18;
}

namespace libMesh
{




// ------------------------------------------------------------
// Tri7 class static member initializations
const int Tri7::num_nodes;
const int Tri7::num_sides;
const int Tri7::num_children;
const int Tri7::nodes_per_side;

const unsigned int Tri7::side_nodes_map[Tri7::num_sides][Tri7::nodes_per_side] =
  {
    {0, 1, 3}, // Side 0
    {1, 2, 4}, // Side 1
    {2, 0, 5}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

const float Tri7::_embedding_matrix[Tri7::num_children][Tri7::num_nodes][Tri7::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //  0      1      2     3      4     5      6
      {   1.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0}, // 0
      {   0.0,   0.0,   0.0,  1.0,   0.0,  0.0,   0.0}, // 1
      {   0.0,   0.0,   0.0,  0.0,   0.0,  1.0,   0.0}, // 2
      {  .375, -.125,   0.0,  .75,   0.0,  0.0,   0.0}, // 3
      {.09375,-.03125,-.03125, .125, -.125, .125,.84375}, // 4
      {  .375,   0.0, -.125,  0.0,   0.0,  .75,   0.0}, // 5
      { 5/r18,-1/r18,-1/r18,4/r18,-2/r18,4/r18,   0.5}  // 6
    },

    // embedding matrix for child 1
    {
      //  0      1      2     3     4      5      6
      {   0.0,   0.0,   0.0,  1.0,  0.0,   0.0,   0.0}, // 0
      {   0.0,   1.0,   0.0,  0.0,  0.0,   0.0,   0.0}, // 1
      {   0.0,   0.0,   0.0,  0.0,  1.0,   0.0,   0.0}, // 2
      { -.125,  .375,   0.0,  .75,  0.0,   0.0,   0.0}, // 3
      {   0.0,  .375, -.125,  0.0,  .75,   0.0,   0.0}, // 4
      {-.03125,.09375,-.03125, .125, .125, -.125,.84375}, // 5
      {-1/r18, 5/r18,-1/r18,4/r18,4/r18,-2/r18,   0.5}  // 6
    },

    // embedding matrix for child 2
    {
      //  0       1     2      3     4     5    6
      {   0.0,   0.0,   0.0,   0.0,  0.0,  1.0,   0.0}, // 0
      {   0.0,   0.0,   0.0,   0.0,  1.0,  0.0,   0.0}, // 1
      {   0.0,   0.0,   1.0,   0.0,  0.0,  0.0,   0.0}, // 2
      {-.03125,-.03125,.09375, -.125, .125, .125,.84375}, // 3
      {   0.0, -.125,  .375,   0.0,  .75,  0.0,   0.0}, // 4
      { -.125,   0.0,  .375,   0.0,  0.0,  .75,   0.0}, // 5
      {-1/r18,-1/r18, 5/r18,-2/r18,4/r18,4/r18,   0.5}  // 6
    },

    // embedding matrix for child 3
    {
      //  0      1      2     3     4     5      6
      {   0.0,   0.0,   0.0,  1.0,  0.0,  0.0,   0.0}, // 0
      {   0.0,   0.0,   0.0,  0.0,  1.0,  0.0,   0.0}, // 1
      {   0.0,   0.0,   0.0,  0.0,  0.0,  1.0,   0.0}, // 2
      {-.03125,.09375,-.03125, .125, .125,-.125,.84375}, // 3
      {-.03125,-.03125,.09375,-.125, .125, .125,.84375}, // 4
      {.09375,-.03125,-.03125, .125,-.125, .125,.84375}, // 5
      {   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   1.0}  // 6
    }
  };

const std::vector<std::pair<unsigned char, unsigned char>>
  Tri7::_parent_bracketing_nodes[Tri7::num_children][Tri7::num_nodes] =
  {
    // Child 0
    {    {{}},   {{}},{{0,2}},{{0,3}},{{3,5}},{{0,5}},{{0,6}} },
    // Child 1
    { {{0,1}},   {{}},{{1,2}},{{1,3}},{{1,4}},{{3,4}},{{1,6}} },
    // Child 2
    { {{0,2}},{{1,2}},   {{}},{{4,5}},{{2,4}},{{2,5}},{{2,6}} },
    // Child 3
    { {{0,1}},{{1,2}},{{0,2}},{{3,4}},{{4,5}},{{3,5}},   {{}} }
  };
#endif



// ------------------------------------------------------------
// Tri7 class member functions

bool Tri7::is_vertex(const unsigned int i) const
{
  if (i < 3)
    return true;
  return false;
}

bool Tri7::is_edge(const unsigned int i) const
{
  if (i < 3)
    return false;
  return true;
}

bool Tri7::is_face(const unsigned int) const
{
  return false;
}

bool Tri7::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Tri7::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Tri7::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
}

bool Tri7::has_affine_map() const
{
  // Make sure edges are straight
  Point v = this->point(2) - this->point(1);
  if (!v.relative_fuzzy_equals
      ((this->point(4) - this->point(1))*2))
    return false;
  v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals
      ((this->point(3) - this->point(0))*2))
    return false;
  Point v20 = this->point(2) - this->point(0);
  if (!v20.relative_fuzzy_equals
      ((this->point(5) - this->point(0))*2))
    return false;

  // Make sure center node is centered
  v += v20;
  if (!v.relative_fuzzy_equals
      ((this->point(6) - this->point(0))*3))
    return false;

  return true;
}



Order Tri7::default_order() const
{
  return THIRD;
}



dof_id_type Tri7::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:

      return
        this->compute_key (this->node_id(3));

    case 1:

      return
        this->compute_key (this->node_id(4));

    case 2:

      return
        this->compute_key (this->node_id(5));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int Tri7::local_side_node(unsigned int side,
                                   unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Tri7::nodes_per_side);

  return Tri7::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> Tri7::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  return this->simple_build_side_ptr<Edge3, Tri7>(i, proxy);
}



void Tri7::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->simple_build_side_ptr<Tri7>(side, i, EDGE3);
}



void Tri7::connectivity(const unsigned int sf,
                        const IOPackage iop,
                        std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        switch(sf)
          {
          case 0:
            // linear sub-triangle 0
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(3)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(5)+1;

            return;

          case 1:
            // linear sub-triangle 1
            conn[0] = this->node_id(3)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(4)+1;
            conn[3] = this->node_id(4)+1;

            return;

          case 2:
            // linear sub-triangle 2
            conn[0] = this->node_id(5)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(2)+1;
            conn[3] = this->node_id(2)+1;

            return;

          case 3:
            // linear sub-triangle 3
            conn[0] = this->node_id(3)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(5)+1;

            return;

          default:
            libmesh_error_msg("Invalid sf = " << sf);
          }
      }

    case VTK:
      {
        // VTK_QUADRATIC_TRIANGLE has same numbering as libmesh TRI6
        conn.resize(6);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        conn[3] = this->node_id(3);
        conn[4] = this->node_id(4);
        conn[5] = this->node_id(5);
        return;

        // Used to write out linear sub-triangles for VTK...
        /*
          conn.resize(3);
          switch(sf)
          {
          case 0:
          // linear sub-triangle 0
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(3);
          conn[2] = this->node_id(5);

          return;

          case 1:
          // linear sub-triangle 1
          conn[0] = this->node_id(3);
          conn[1] = this->node_id(1);
          conn[2] = this->node_id(4);

          return;

          case 2:
          // linear sub-triangle 2
          conn[0] = this->node_id(5);
          conn[1] = this->node_id(4);
          conn[2] = this->node_id(2);

          return;

          case 3:
          // linear sub-triangle 3
          conn[0] = this->node_id(3);
          conn[1] = this->node_id(4);
          conn[2] = this->node_id(5);

          return;

          default:
          libmesh_error_msg("Invalid sf = " << sf);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



BoundingBox Tri7::loose_bounding_box () const
{
  // This might have curved edges, or might be a curved surface in
  // 3-space, in which case the full bounding box can be larger than
  // the bounding box of just the nodes.
  //
  //
  // FIXME - I haven't yet proven the formula below to be correct for
  // quadratics in 2D - RHS
  //
  // FIXME - This doesn't take into account curvature caused by the
  // center node in 3D - RHS
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      Real center = this->point(0)(d);
      for (unsigned int p=1; p != 6; ++p)
        center += this->point(p)(d);
      center /= 6;

      Real hd = std::abs(center - this->point(0)(d));
      for (unsigned int p=1; p != 6; ++p)
        hd = std::max(hd, std::abs(center - this->point(p)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



unsigned int Tri7::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 3:
    case 4:
    case 5:
      return 2;

    case 6:
      return 3;

    default:
      libmesh_error_msg("Invalid n = " << n);
    }
}



unsigned short int Tri7::second_order_adjacent_vertex (const unsigned int n,
                                                       const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
    case 6:
      {
        libmesh_assert_less (v, 3);
        return static_cast<unsigned short int>(v);
      }

    default:
      {
        libmesh_assert_less (v, 2);
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }
    }
}



const unsigned short int Tri7::_second_order_adjacent_vertices[Tri7::num_sides][2] =
  {
    {0, 1}, // vertices adjacent to node 3
    {1, 2}, // vertices adjacent to node 4
    {0, 2}  // vertices adjacent to node 5
  };



std::pair<unsigned short int, unsigned short int>
Tri7::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



const unsigned short int Tri7::_second_order_vertex_child_number[Tri7::num_nodes] =
  {
    99,99,99, // Vertices
    0,1,0,    // Edges
    3         // Interior
  };



const unsigned short int Tri7::_second_order_vertex_child_index[Tri7::num_nodes] =
  {
    99,99,99, // Vertices
    1,2,2,    // Edges
    6         // Interior
  };


void Tri7::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 3);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
    }
}


unsigned int Tri7::center_node_on_side(const unsigned short side) const
{
  libmesh_assert_less (side, Tri7::num_sides);
  return side + 3;
}


} // namespace libMesh
