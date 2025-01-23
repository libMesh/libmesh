// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_prism21.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_tri7.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"

namespace {
  // Avoid downcasting when Real > double
  constexpr libMesh::Real r6 = 6;
  constexpr libMesh::Real r9 = 9;
  constexpr libMesh::Real r12 = 12;
  constexpr libMesh::Real r18 = 18;
  constexpr libMesh::Real r24 = 24;
  constexpr libMesh::Real r36 = 36;
  constexpr libMesh::Real r48 = 48;
  constexpr libMesh::Real r72 = 72;
  constexpr libMesh::Real r144 = 144;
}

namespace libMesh
{



// ------------------------------------------------------------
// Prism21 class static member initializations
const int Prism21::num_nodes;
const int Prism21::nodes_per_side;
const int Prism21::nodes_per_edge;

const unsigned int Prism21::side_nodes_map[Prism21::num_sides][Prism21::nodes_per_side] =
  {
    {0, 2, 1,  8,  7,  6, 18, 99, 99}, // Side 0
    {0, 1, 4,  3,  6, 10, 12,  9, 15}, // Side 1
    {1, 2, 5,  4,  7, 11, 13, 10, 16}, // Side 2
    {2, 0, 3,  5,  8,  9, 14, 11, 17}, // Side 3
    {3, 4, 5, 12, 13, 14, 19, 99, 99}  // Side 4
  };

const unsigned int Prism21::edge_nodes_map[Prism21::num_edges][Prism21::nodes_per_edge] =
  {
    {0, 1,  6}, // Edge 0
    {1, 2,  7}, // Edge 1
    {0, 2,  8}, // Edge 2
    {0, 3,  9}, // Edge 3
    {1, 4, 10}, // Edge 4
    {2, 5, 11}, // Edge 5
    {3, 4, 12}, // Edge 6
    {4, 5, 13}, // Edge 7
    {3, 5, 14}  // Edge 8
  };

// ------------------------------------------------------------
// Prism21 class member functions

bool Prism21::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

bool Prism21::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  if (i > 14)
    return false;
  return true;
}

bool Prism21::is_face(const unsigned int i) const
{
  if (i > 19)
    return false;
  if (i > 14)
    return true;
  return false;
}

bool Prism21::is_node_on_side(const unsigned int n,
                              const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Prism21::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s > 0 && s < 4) ? 0 : 2;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
Prism21::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Prism21::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



bool Prism21::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2), affine_tol))
    return false;

  // Make sure edges are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(9) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(11) - this->point(2), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(15) - this->point(6), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(16) - this->point(7), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(17) - this->point(8), affine_tol))
    return false;
  v = (this->point(1) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(6) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(12) - this->point(3), affine_tol))
    return false;
  v = (this->point(2) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(3), affine_tol))
    return false;
  v = (this->point(2) - this->point(1))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(4), affine_tol))
    return false;

  // Make sure triangle face midpoints are centered
  v = this->point(2) + this->point(1) - 2*this->point(0);
  if (!v.relative_fuzzy_equals((this->point(18)-this->point(0))*3))
    return false;
  v = this->point(5) + this->point(4) - 2*this->point(3);
  if (!v.relative_fuzzy_equals((this->point(19)-this->point(3))*3))
    return false;
  v = this->point(11) + this->point(10) - 2*this->point(9);
  if (!v.relative_fuzzy_equals((this->point(20)-this->point(9))*3))
    return false;

  return true;
}



Order Prism21::default_order() const
{
  return THIRD;
}

dof_id_type Prism21::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:  // the triangular face at z=0
      {
        return Elem::compute_key (this->node_id(18));
      }
    case 1:  // the quad face at y=0
      {
        return Elem::compute_key (this->node_id(15));
      }
    case 2:  // the other quad face
      {
        return Elem::compute_key (this->node_id(16));
      }
    case 3: // the quad face at x=0
      {
        return Elem::compute_key (this->node_id(17));
      }
    case 4: // the triangular face at z=1
      {
        return Elem::compute_key (this->node_id(19));
      }
    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



unsigned int Prism21::local_side_node(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 9 nodes per side.
  libmesh_assert_less(side_node, Prism21::nodes_per_side);

  // Some sides have 7 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 7);

  return Prism21::side_nodes_map[side][side_node];
}



unsigned int Prism21::local_edge_node(unsigned int edge,
                                      unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, Prism21::nodes_per_edge);

  return Prism21::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Prism21::build_side_ptr (const unsigned int i,
                                               bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;
  if (proxy)
    {
#ifdef LIBMESH_ENABLE_DEPRECATED
      libmesh_deprecated();
      switch(i)
        {
        case 0:
        case 4:
          {
            face = std::make_unique<Side<Tri7,Prism21>>(this,i);
            break;
          }

        case 1:
        case 2:
        case 3:
          {
            face = std::make_unique<Side<Quad9,Prism21>>(this,i);
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
#else
      libmesh_error();
#endif // LIBMESH_ENABLE_DEPRECATED
    }
  else
    {
      switch (i)
        {
        case 0: // the triangular face at z=-1
        case 4: // the triangular face at z=1
          {
            face = std::make_unique<Tri7>();
            break;
          }
        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = std::make_unique<Quad9>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Prism21::side_nodes_map[i][n]);
    }

#ifdef LIBMESH_ENABLE_DEPRECATED
  if (!proxy) // proxy sides used to leave parent() set
#endif
    face->set_parent(nullptr);
  face->set_interior_parent(this);

  face->subdomain_id() = this->subdomain_id();
  face->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  face->set_p_level(this->p_level());
#endif

  return face;
}



void Prism21::build_side_ptr (std::unique_ptr<Elem> & side,
                              const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // the triangular face at z=-1
    case 4: // the triangular face at z=1
      {
        if (!side.get() || side->type() != TRI7)
          {
            side = this->build_side_ptr(i, false);
            return;
          }
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        if (!side.get() || side->type() != QUAD9)
          {
            side = this->build_side_ptr(i, false);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->subdomain_id() = this->subdomain_id();
  side->set_mapping_type(this->mapping_type());

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n) = this->node_ptr(Prism21::side_nodes_map[i][n]);
}



std::unique_ptr<Elem> Prism21::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Prism21>(i);
}



void Prism21::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Prism21>(edge, i, EDGE3);
}



void Prism21::connectivity(const unsigned int /*sc*/,
                           const IOPackage /*iop*/,
                           std::vector<dof_id_type> & /*conn*/) const
{
  libmesh_not_implemented(); // FIXME RHS

/*
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        switch (sc)
          {

          case 0:
            {
              conn[0] = this->node_id(0)+1;
              conn[1] = this->node_id(6)+1;
              conn[2] = this->node_id(8)+1;
              conn[3] = this->node_id(8)+1;
              conn[4] = this->node_id(9)+1;
              conn[5] = this->node_id(15)+1;
              conn[6] = this->node_id(17)+1;
              conn[7] = this->node_id(17)+1;

              return;
            }

          case 1:
            {
              conn[0] = this->node_id(6)+1;
              conn[1] = this->node_id(1)+1;
              conn[2] = this->node_id(7)+1;
              conn[3] = this->node_id(7)+1;
              conn[4] = this->node_id(15)+1;
              conn[5] = this->node_id(10)+1;
              conn[6] = this->node_id(16)+1;
              conn[7] = this->node_id(16)+1;

              return;
            }

          case 2:
            {
              conn[0] = this->node_id(8)+1;
              conn[1] = this->node_id(7)+1;
              conn[2] = this->node_id(2)+1;
              conn[3] = this->node_id(2)+1;
              conn[4] = this->node_id(17)+1;
              conn[5] = this->node_id(16)+1;
              conn[6] = this->node_id(11)+1;
              conn[7] = this->node_id(11)+1;

              return;
            }

          case 3:
            {
              conn[0] = this->node_id(6)+1;
              conn[1] = this->node_id(7)+1;
              conn[2] = this->node_id(8)+1;
              conn[3] = this->node_id(8)+1;
              conn[4] = this->node_id(15)+1;
              conn[5] = this->node_id(16)+1;
              conn[6] = this->node_id(17)+1;
              conn[7] = this->node_id(17)+1;

              return;
            }

          case 4:
            {
              conn[0] = this->node_id(9)+1;
              conn[1] = this->node_id(15)+1;
              conn[2] = this->node_id(17)+1;
              conn[3] = this->node_id(17)+1;
              conn[4] = this->node_id(3)+1;
              conn[5] = this->node_id(12)+1;
              conn[6] = this->node_id(14)+1;
              conn[7] = this->node_id(14)+1;

              return;
            }

          case 5:
            {
              conn[0] = this->node_id(15)+1;
              conn[1] = this->node_id(10)+1;
              conn[2] = this->node_id(16)+1;
              conn[3] = this->node_id(16)+1;
              conn[4] = this->node_id(12)+1;
              conn[5] = this->node_id(4)+1;
              conn[6] = this->node_id(13)+1;
              conn[7] = this->node_id(13)+1;

              return;
            }

          case 6:
            {
              conn[0] = this->node_id(17)+1;
              conn[1] = this->node_id(16)+1;
              conn[2] = this->node_id(11)+1;
              conn[3] = this->node_id(11)+1;
              conn[4] = this->node_id(14)+1;
              conn[5] = this->node_id(13)+1;
              conn[6] = this->node_id(5)+1;
              conn[7] = this->node_id(5)+1;

              return;
            }

          case 7:
            {
              conn[0] = this->node_id(15)+1;
              conn[1] = this->node_id(16)+1;
              conn[2] = this->node_id(17)+1;
              conn[3] = this->node_id(17)+1;
              conn[4] = this->node_id(12)+1;
              conn[5] = this->node_id(13)+1;
              conn[6] = this->node_id(14)+1;
              conn[7] = this->node_id(14)+1;

              return;
            }

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }

      }

    case VTK:
      {
        // VTK now supports VTK_BIQUADRATIC_QUADRATIC_WEDGE directly
        const unsigned int conn_size = 18;
        conn.resize(conn_size);

        // VTK's VTK_BIQUADRATIC_QUADRATIC_WEDGE first 9 (vertex) and
        // last 3 (mid-face) nodes match.  The middle and top layers
        // of mid-edge nodes are reversed from LibMesh's.
        for (auto i : index_range(conn))
          conn[i] = this->node_id(i);

        // top "ring" of mid-edge nodes
        conn[9]  = this->node_id(12);
        conn[10] = this->node_id(13);
        conn[11] = this->node_id(14);

        // middle "ring" of mid-edge nodes
        conn[12] = this->node_id(9);
        conn[13] = this->node_id(10);
        conn[14] = this->node_id(11);

        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
*/
}




unsigned int Prism21::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
      return 2;

    case 15:
    case 16:
    case 17:
      return 4;

    case 18:
    case 19:
      return 3;

    case 20:
      return 6;

    default:
      libmesh_error_msg("Invalid node n = " << n);
    }
}





unsigned short int Prism21::second_order_adjacent_vertex (const unsigned int n,
                                                          const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
      /*
       * These nodes are unique to \p Prism20,
       * let our _remaining_... matrix handle
       * this.
       */
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:
      {
        libmesh_assert_less (v, 4);
        return _remaining_second_order_adjacent_vertices[n-15][v];
      }

      /*
       * for the bubble node the return value is simply v.
       * Why? -- the user asks for the v-th adjacent vertex,
       * from \p n_second_order_adjacent_vertices() there
       * are 6 adjacent vertices, and these happen to be
       * 0..5
       */
    case 20:
      {
        libmesh_assert_less (v, 6);
        return static_cast<unsigned short int>(v);
      }

      /*
       * All other second-order nodes (6,...,14) are
       * identical with Prism15 and are therefore
       * delegated to the _second_order matrix of
       * \p Prism
       */
    default:
      {
        libmesh_assert_less (v, 2);
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }

    }

  return static_cast<unsigned short int>(-1);
}



const unsigned short int Prism21::_remaining_second_order_adjacent_vertices[5][4] =
  {
    { 0,  1,  3,  4}, // vertices adjacent to node 15
    { 1,  2,  4,  5}, // vertices adjacent to node 16
    { 0,  2,  3,  5}, // vertices adjacent to node 17
    { 0,  1,  2, 99}, // vertices adjacent to node 18
    { 3,  4,  5, 99}  // vertices adjacent to node 19
  };



std::pair<unsigned short int, unsigned short int>
Prism21::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



#ifdef LIBMESH_ENABLE_AMR

// FIXME RHS

const Real Prism21::_embedding_matrix[Prism21::num_children][Prism21::num_nodes][Prism21::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 5
      {           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 6
      {          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 7
      {           3/8.,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 8
      {           3/8.,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0}, // 11
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 12
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,         27/32.}, // 13
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0}, // 14
      {          9/64.,         -3/64.,              0,         -3/64.,          1/64.,              0,          9/32.,              0,              0,          9/32.,         -3/32.,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0,              0,              0}, // 15
      {         9/256.,        -3/256.,        -3/256.,        -3/256.,         1/256.,         1/256.,          3/64.,         -3/64.,          3/64.,         9/128.,        -3/128.,        -3/128.,         -1/64.,          1/64.,         -1/64.,          3/32.,         -3/32.,          3/32.,        81/256.,       -27/256.,        81/128.}, // 16
      {          9/64.,              0,         -3/64.,         -3/64.,              0,          1/64.,              0,              0,          9/32.,          9/32.,              0,         -3/32.,              0,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0}, // 17
      {          5/r18,         -1/r18,         -1/r18,              0,              0,              0,           2/r9,          -1/r9,           2/r9,              0,              0,              0,              0,              0,              0,              0,              0,              0,           1/2.,              0,              0}, // 18
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          5/r18,         -1/r18,         -1/r18,              0,              0,              0,           2/r9,          -1/r9,           2/r9,              0,              0,           1/2.}, // 19
      {          5/r48,         -1/r48,         -1/r48,        -5/r144,         1/r144,         1/r144,          1/r12,         -1/r24,          1/r12,          5/r24,         -1/r24,         -1/r24,         -1/r36,          1/r72,         -1/r36,           1/r6,         -1/r12,           1/r6,          3/16.,         -1/16.,           3/8.}  // 20
    },
    // embedding matrix for child 1
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 0
      {              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 5
      {          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 6
      {              0,           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 7
      {         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 8
      {              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 9
      {              0,           3/8.,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 11
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 12
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0}, // 13
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,         27/32.}, // 14
      {         -3/64.,          9/64.,              0,          1/64.,         -3/64.,              0,          9/32.,              0,              0,         -3/32.,          9/32.,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0,              0,              0}, // 15
      {              0,          9/64.,         -3/64.,              0,         -3/64.,          1/64.,              0,          9/32.,              0,              0,          9/32.,         -3/32.,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0,              0}, // 16
      {        -3/256.,         9/256.,        -3/256.,         1/256.,        -3/256.,         1/256.,          3/64.,          3/64.,         -3/64.,        -3/128.,         9/128.,        -3/128.,         -1/64.,         -1/64.,          1/64.,          3/32.,          3/32.,         -3/32.,        81/256.,       -27/256.,        81/128.}, // 17
      {         -1/r18,          5/r18,         -1/r18,              0,              0,              0,           2/r9,           2/r9,          -1/r9,              0,              0,              0,              0,              0,              0,              0,              0,              0,           1/2.,              0,              0}, // 18
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/r18,          5/r18,         -1/r18,              0,              0,              0,           2/r9,           2/r9,          -1/r9,              0,              0,           1/2.}, // 19
      {         -1/r48,          5/r48,         -1/r48,         1/r144,        -5/r144,         1/r144,          1/r12,          1/r12,         -1/r24,         -1/r24,          5/r24,         -1/r24,         -1/r36,         -1/r36,          1/r72,           1/r6,           1/r6,         -1/r12,          3/16.,         -1/16.,           3/8.}  // 20
    },
    // embedding matrix for child 2
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 5
      {         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 6
      {              0,          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 7
      {          -1/8.,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 8
      {              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 10
      {              0,              0,           3/8.,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 11
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,         27/32.}, // 12
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0}, // 13
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0}, // 14
      {        -3/256.,        -3/256.,         9/256.,         1/256.,         1/256.,        -3/256.,         -3/64.,          3/64.,          3/64.,        -3/128.,        -3/128.,         9/128.,          1/64.,         -1/64.,         -1/64.,         -3/32.,          3/32.,          3/32.,        81/256.,       -27/256.,        81/128.}, // 15
      {              0,         -3/64.,          9/64.,              0,          1/64.,         -3/64.,              0,          9/32.,              0,              0,         -3/32.,          9/32.,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0,              0}, // 16
      {         -3/64.,              0,          9/64.,          1/64.,              0,         -3/64.,              0,              0,          9/32.,         -3/32.,              0,          9/32.,              0,              0,         -3/32.,              0,              0,          9/16.,              0,              0,              0}, // 17
      {         -1/r18,         -1/r18,          5/r18,              0,              0,              0,          -1/r9,           2/r9,           2/r9,              0,              0,              0,              0,              0,              0,              0,              0,              0,           1/2.,              0,              0}, // 18
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/r18,         -1/r18,          5/r18,              0,              0,              0,          -1/r9,           2/r9,           2/r9,              0,              0,           1/2.}, // 19
      {         -1/r48,         -1/r48,          5/r48,         1/r144,         1/r144,        -5/r144,         -1/r24,          1/r12,          1/r12,         -1/r24,         -1/r24,          5/r24,          1/r72,         -1/r36,         -1/r36,         -1/r12,           1/r6,           1/r6,          3/16.,         -1/16.,           3/8.}  // 20
    },
    // embedding matrix for child 3
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 5
      {         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 6
      {         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 7
      {          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,              0,              0,              0,              0,              0,              0,              0,         27/32.,              0,              0}, // 8
      {              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,              0,              0,              0,              0,          -1/8.,              0,              0,           3/4.,              0,              0,              0}, // 11
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,         27/32.}, // 12
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,         27/32.}, // 13
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,         27/32.}, // 14
      {        -3/256.,         9/256.,        -3/256.,         1/256.,        -3/256.,         1/256.,          3/64.,          3/64.,         -3/64.,        -3/128.,         9/128.,        -3/128.,         -1/64.,         -1/64.,          1/64.,          3/32.,          3/32.,         -3/32.,        81/256.,       -27/256.,        81/128.}, // 15
      {        -3/256.,        -3/256.,         9/256.,         1/256.,         1/256.,        -3/256.,         -3/64.,          3/64.,          3/64.,        -3/128.,        -3/128.,         9/128.,          1/64.,         -1/64.,         -1/64.,         -3/32.,          3/32.,          3/32.,        81/256.,       -27/256.,        81/128.}, // 16
      {         9/256.,        -3/256.,        -3/256.,        -3/256.,         1/256.,         1/256.,          3/64.,         -3/64.,          3/64.,         9/128.,        -3/128.,        -3/128.,         -1/64.,          1/64.,         -1/64.,          3/32.,         -3/32.,          3/32.,        81/256.,       -27/256.,        81/128.}, // 17
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0}, // 18
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1}, // 19
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,          -1/8.,           3/4.}  // 20
    },
    // embedding matrix for child 4
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 2
      {              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0}, // 5
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 6
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,         27/32.}, // 7
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,              0,          -1/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0}, // 8
      {          -1/8.,              0,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0}, // 11
      {              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0}, // 12
      {              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,              0,              0,         27/32.,              0}, // 13
      {              0,              0,              0,           3/8.,              0,          -1/8.,              0,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0}, // 14
      {         -3/64.,          1/64.,              0,          9/64.,         -3/64.,              0,         -3/32.,              0,              0,          9/32.,         -3/32.,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0,              0,              0}, // 15
      {        -3/256.,         1/256.,         1/256.,         9/256.,        -3/256.,        -3/256.,         -1/64.,          1/64.,         -1/64.,         9/128.,        -3/128.,        -3/128.,          3/64.,         -3/64.,          3/64.,          3/32.,         -3/32.,          3/32.,       -27/256.,        81/256.,        81/128.}, // 16
      {         -3/64.,              0,          1/64.,          9/64.,              0,         -3/64.,              0,              0,         -3/32.,          9/32.,              0,         -3/32.,              0,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0}, // 17
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          5/r18,         -1/r18,         -1/r18,              0,              0,              0,           2/r9,          -1/r9,           2/r9,              0,              0,           1/2.}, // 18
      {              0,              0,              0,          5/r18,         -1/r18,         -1/r18,              0,              0,              0,              0,              0,              0,           2/r9,          -1/r9,           2/r9,              0,              0,              0,              0,           1/2.,              0}, // 19
      {        -5/r144,         1/r144,         1/r144,          5/r48,         -1/r48,         -1/r48,         -1/r36,          1/r72,         -1/r36,          5/r24,         -1/r24,         -1/r24,          1/r12,         -1/r24,          1/r12,           1/r6,         -1/r12,           1/r6,         -1/16.,          3/16.,           3/8.}  // 20
    },
    // embedding matrix for child 5
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0}, // 5
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 6
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0}, // 7
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,         27/32.}, // 8
      {              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 9
      {              0,          -1/8.,              0,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 11
      {              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0}, // 12
      {              0,              0,              0,              0,           3/8.,          -1/8.,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0}, // 13
      {              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,              0,              0,         27/32.,              0}, // 14
      {          1/64.,         -3/64.,              0,         -3/64.,          9/64.,              0,         -3/32.,              0,              0,         -3/32.,          9/32.,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0,              0,              0}, // 15
      {              0,         -3/64.,          1/64.,              0,          9/64.,         -3/64.,              0,         -3/32.,              0,              0,          9/32.,         -3/32.,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0,              0}, // 16
      {         1/256.,        -3/256.,         1/256.,        -3/256.,         9/256.,        -3/256.,         -1/64.,         -1/64.,          1/64.,        -3/128.,         9/128.,        -3/128.,          3/64.,          3/64.,         -3/64.,          3/32.,          3/32.,         -3/32.,       -27/256.,        81/256.,        81/128.}, // 17
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/r18,          5/r18,         -1/r18,              0,              0,              0,           2/r9,           2/r9,          -1/r9,              0,              0,           1/2.}, // 18
      {              0,              0,              0,         -1/r18,          5/r18,         -1/r18,              0,              0,              0,              0,              0,              0,           2/r9,           2/r9,          -1/r9,              0,              0,              0,              0,           1/2.,              0}, // 19
      {         1/r144,        -5/r144,         1/r144,         -1/r48,          5/r48,         -1/r48,         -1/r36,         -1/r36,          1/r72,         -1/r24,          5/r24,         -1/r24,          1/r12,          1/r12,         -1/r24,           1/r6,           1/r6,         -1/r12,         -1/16.,          3/16.,           3/8.}  // 20
    },
    // embedding matrix for child 6
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 5
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,         27/32.}, // 6
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0}, // 7
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0}, // 8
      {              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 10
      {              0,              0,          -1/8.,              0,              0,           3/8.,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0,              0,              0}, // 11
      {              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,              0,              0,         27/32.,              0}, // 12
      {              0,              0,              0,              0,          -1/8.,           3/8.,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0,              0}, // 13
      {              0,              0,              0,          -1/8.,              0,           3/8.,              0,              0,              0,              0,              0,              0,              0,              0,           3/4.,              0,              0,              0,              0,              0,              0}, // 14
      {         1/256.,         1/256.,        -3/256.,        -3/256.,        -3/256.,         9/256.,          1/64.,         -1/64.,         -1/64.,        -3/128.,        -3/128.,         9/128.,         -3/64.,          3/64.,          3/64.,         -3/32.,          3/32.,          3/32.,       -27/256.,        81/256.,        81/128.}, // 15
      {              0,          1/64.,         -3/64.,              0,         -3/64.,          9/64.,              0,         -3/32.,              0,              0,         -3/32.,          9/32.,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0,              0}, // 16
      {          1/64.,              0,         -3/64.,         -3/64.,              0,          9/64.,              0,              0,         -3/32.,         -3/32.,              0,          9/32.,              0,              0,          9/32.,              0,              0,          9/16.,              0,              0,              0}, // 17
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/r18,         -1/r18,          5/r18,              0,              0,              0,          -1/r9,           2/r9,           2/r9,              0,              0,           1/2.}, // 18
      {              0,              0,              0,         -1/r18,         -1/r18,          5/r18,              0,              0,              0,              0,              0,              0,          -1/r9,           2/r9,           2/r9,              0,              0,              0,              0,           1/2.,              0}, // 19
      {         1/r144,         1/r144,        -5/r144,         -1/r48,         -1/r48,          5/r48,          1/r72,         -1/r36,         -1/r36,         -1/r24,         -1/r24,          5/r24,         -1/r24,          1/r12,          1/r12,         -1/r12,           1/r6,           1/r6,         -1/16.,          3/16.,           3/8.}  // 20
    },
    // embedding matrix for child 7
    {
      //             0,               1,              2,              3,              4,              5,              6,              7,              8,              9,              10,             11,             12,             13,             14,             15,             16,             17,             18,             19,             20
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0}, // 0
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0}, // 1
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0}, // 2
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0,              0}, // 3
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0,              0}, // 4
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0,              0,              0,              0,              0,              0}, // 5
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,         27/32.}, // 6
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,         27/32.}, // 7
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,         27/32.}, // 8
      {              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0,              0}, // 9
      {              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0,              0}, // 10
      {              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,              0,              0,              0,              0,              0,           3/8.,              0,              0,           3/4.,              0,              0,              0}, // 11
      {              0,              0,              0,         -1/32.,          3/32.,         -1/32.,              0,              0,              0,              0,              0,              0,           1/8.,           1/8.,          -1/8.,              0,              0,              0,              0,         27/32.,              0}, // 12
      {              0,              0,              0,         -1/32.,         -1/32.,          3/32.,              0,              0,              0,              0,              0,              0,          -1/8.,           1/8.,           1/8.,              0,              0,              0,              0,         27/32.,              0}, // 13
      {              0,              0,              0,          3/32.,         -1/32.,         -1/32.,              0,              0,              0,              0,              0,              0,           1/8.,          -1/8.,           1/8.,              0,              0,              0,              0,         27/32.,              0}, // 14
      {         1/256.,        -3/256.,         1/256.,        -3/256.,         9/256.,        -3/256.,         -1/64.,         -1/64.,          1/64.,        -3/128.,         9/128.,        -3/128.,          3/64.,          3/64.,         -3/64.,          3/32.,          3/32.,         -3/32.,       -27/256.,        81/256.,        81/128.}, // 15
      {         1/256.,         1/256.,        -3/256.,        -3/256.,        -3/256.,         9/256.,          1/64.,         -1/64.,         -1/64.,        -3/128.,        -3/128.,         9/128.,         -3/64.,          3/64.,          3/64.,         -3/32.,          3/32.,          3/32.,       -27/256.,        81/256.,        81/128.}, // 16
      {        -3/256.,         1/256.,         1/256.,         9/256.,        -3/256.,        -3/256.,         -1/64.,          1/64.,         -1/64.,         9/128.,        -3/128.,        -3/128.,          3/64.,         -3/64.,          3/64.,          3/32.,         -3/32.,          3/32.,       -27/256.,        81/256.,        81/128.}, // 17
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1}, // 18
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              1,              0}, // 19
      {              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,          -1/8.,           3/8.,           3/4.}  // 20
    }
  };

const std::vector<std::pair<unsigned char, unsigned char>>
  Prism21::_parent_bracketing_nodes[Prism21::num_children][Prism21::num_nodes] =
  {
    // Child 0
    {      {},{{0,1}},{{0,2}},{{0,3}},{{0,4},{1,3},{6,12},{9,10}},{{0,5},{2,3},{8,14},{9,11}},
      {{0,6}},{{6,8}},{{0,8}},{{0,9}},{{6,15}},{{8,17}},{{9,15}},{{15,17}},{{9,17}},
      {{0,15},{6,9}},{{6,17},{8,15}},{{0,17},{8,9}},{{0,18}},{{9,20}},{{0,20},{9,18}} },
    // Child 1
    { {{0,1}},     {},{{1,2}},{{0,4},{1,3},{6,12},{9,10}},{{1,4}},{{1,5},{2,4},{7,13},{10,11}},
      {{1,6}},{{1,7}},{{6,7}},{{6,15}},{{1,10}},{{7,16}},{{10,15}},{{10,16}},{{15,16}},
      {{1,15},{6,10}},{{1,16},{7,10}},{{6,16},{7,15}},{{1,18}},{{10,20}},{{1,20},{10,18}} },
    // Child 2
    { {{0,2}},{{1,2}},     {},{{0,5},{2,3},{8,14},{9,11}},{{1,5},{2,4},{7,13},{10,11}},{{2,5}},
      {{7,8}},{{2,7}},{{2,8}},{{8,17}},{{7,16}},{{2,11}},{{16,17}},{{11,16}},{{11,17}},
      {{7,17},{8,16}},{{2,16},{7,11}},{{2,17},{8,11}},{{2,18}},{{11,20}},{{2,20},{11,18}} },
    // Child 3
    { {{0,1}},{{1,2}},{{0,2}},{{0,4},{1,3},{6,12},{9,10}},{{1,5},{2,4},{7,13},{10,11}},{{0,5},{2,3},{8,14},{9,11}},
      {{6,7}},{{7,8}},{{6,8}},{{6,15}},{{7,16}},{{8,17}},{{15,16}},{{16,17}},{{15,17}},
      {{6,16},{7,15}},{{7,17},{8,16}},{{6,17},{8,15}},{{0,7},{1,8},{2,6}},{{9,16},{10,17},{11,15}},{{0,16},{1,17},{2,15}} },
    // Child 4
    { {{0,3}},{{0,4},{1,3},{6,11},{9,10}},{{0,5},{2,3},{8,14},{9,11}},  {},{{3,4}},{{3,5}},
      {{9,15}},{{15,17}},{{9,17}},{{3,9}},{{12,15}},{{14,17}},{{3,12}},{{12,14}},{{3,14}},
      {{3,15},{9,12}},{{12,17},{14,15}},{{3,17},{9,14}},{{9,20}},{{3,19}},{{9,19},{3,20}} },
    // Child 5
    { {{0,4},{1,3},{6,12},{9,10}},{{1,4}},{{1,5},{2,4},{7,13},{10,11}},{{3,4}},  {},{{4,5}},
      {{10,15}},{{10,16}},{{15,16}},{{12,15}},{{4,10}},{{13,16}},{{4,12}},{{4,13}},{{12,13}},
      {{4,15},{10,12}},{{4,16},{10,13}},{{12,16},{13,15}},{{10,20}},{{4,19}},{{10,19},{4,20}} },
    // Child 6
    { {{0,5},{2,3},{8,14},{9,11}},{{1,5},{2,4},{7,13},{10,11}},{{2,5}},{{3,5}},{{4,5}},  {},
      {{16,17}},{{11,16}},{{11,17}},{{14,17}},{{13,16}},{{5,11}},{{13,14}},{{5,13}},{{5,14}},
      {{13,17},{14,16}},{{5,16},{11,13}},{{5,17},{11,14}},{{11,20}},{{5,19}},{{11,19},{5,20}} },
    // Child 7
    { {{0,4},{1,3},{6,12},{9,10}},{{1,5},{2,4},{7,13},{10,11}},{{0,5},{2,3},{8,14},{9,11}},{{3,4}},{{4,5}},{{3,5}},
      {{15,16}},{{16,17}},{{15,17}},{{12,15}},{{13,16}},{{14,17}},{{12,13}},{{13,14}},{{12,14}},
      {{12,16},{13,15}},{{13,17},{14,16}},{{12,17},{14,15}},{{9,16},{10,17},{11,15}},{{3,13},{4,14},{5,12}},{{3,16},{4,17},{5,15}} }
  };

#endif


void
Prism21::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 6);
  const unsigned int side = perm_num % 2;
  const unsigned int rotate = perm_num / 2;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
      swap3nodes(6,7,8);
      swap3nodes(9,10,11);
      swap3nodes(12,13,14);
      swap3nodes(15,16,17);
      swap3neighbors(1,2,3);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap2nodes(1,3);
    swap2nodes(0,4);
    swap2nodes(2,5);
    swap2nodes(6,12);
    swap2nodes(9,10);
    swap2nodes(7,14);
    swap2nodes(8,13);
    swap2nodes(16,17);
    swap2nodes(18,19);
    swap2neighbors(0,4);
    swap2neighbors(2,3);
    break;
  default:
    libmesh_error();
  }
}


void
Prism21::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(3,4);
  swap2nodes(7,8);
  swap2nodes(9,10);
  swap2nodes(13,14);
  swap2nodes(16,17);
  swap2neighbors(2,3);
  swap2boundarysides(2,3,boundary_info);
  swap2boundaryedges(0,1,boundary_info);
  swap2boundaryedges(3,4,boundary_info);
  swap2boundaryedges(7,8,boundary_info);
}


unsigned int Prism21::center_node_on_side(const unsigned short side) const
{
  libmesh_assert_less (side, Prism21::num_sides);
  if (side >= 1 && side <= 3)
    return side + 14;
  if (side == 4)
    return 19;
  return 18;
}


ElemType
Prism21::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s == 0 || s == 4)
    return TRI7;
  return QUAD9;
}


} // namespace libMesh
