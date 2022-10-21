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
#include "libmesh/side.h"
#include "libmesh/cell_prism21.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_tri7.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"

namespace libMesh
{



// ------------------------------------------------------------
// Prism21 class static member initializations
const int Prism21::num_nodes;
const int Prism21::num_sides;
const int Prism21::num_edges;
const int Prism21::num_children;
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



void Prism21::connectivity(const unsigned int sc,
                           const IOPackage iop,
                           std::vector<dof_id_type> & conn) const
{
  libmesh_not_implemented(); // FIXME RHS

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
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
      {    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
      { 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}, // 16
      { 0.140625,       0.,-0.046875,-0.046875,       0., 0.015625,       0.,       0.,  0.28125,  0.28125,       0., -0.09375,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 1
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 5
      {   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 14
      {-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0., 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
      {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}  // 17
    },

    // embedding matrix for child 2
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
      {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
      {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 15
      {       0.,-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
      {-0.046875,       0., 0.140625, 0.015625,       0.,-0.046875,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 3
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
      {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 14
      {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}, // 15
      {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 16
      {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}  // 17
    },

    // embedding matrix for child 4
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
      {       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
      {   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
      {-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}, // 16
      {-0.046875,       0., 0.015625, 0.140625,       0.,-0.046875,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 5
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 11
      {       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 14
      { 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0.,-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
      { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}  // 17
    },

    // embedding matrix for child 6
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
      {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
      { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 15
      {       0., 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
      { 0.015625,       0.,-0.046875,-0.046875,       0., 0.140625,       0.,       0., -0.09375, -0.09375,       0.,  0.28125,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 7
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 14
      { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}, // 15
      { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 16
      {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}  // 17
    }
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
