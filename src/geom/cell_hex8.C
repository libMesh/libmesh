// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"

namespace libMesh
{




// ------------------------------------------------------------
// Hex8 class static member initializations
const unsigned int Hex8::side_nodes_map[6][4] =
  {
    {0, 3, 2, 1}, // Side 0
    {0, 1, 5, 4}, // Side 1
    {1, 2, 6, 5}, // Side 2
    {2, 3, 7, 6}, // Side 3
    {3, 0, 4, 7}, // Side 4
    {4, 5, 6, 7}  // Side 5
  };

const unsigned int Hex8::edge_nodes_map[12][2] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 5}, // Edge 5
    {2, 6}, // Edge 6
    {3, 7}, // Edge 7
    {4, 5}, // Edge 8
    {5, 6}, // Edge 9
    {6, 7}, // Edge 10
    {4, 7}  // Edge 11
  };


// ------------------------------------------------------------
// Hex8 class member functions

bool Hex8::is_vertex(const unsigned int) const
{
  return true;
}

bool Hex8::is_edge(const unsigned int) const
{
  return false;
}

bool Hex8::is_face(const unsigned int) const
{
  return false;
}

bool Hex8::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Hex8::is_node_on_edge(const unsigned int n,
                           const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Hex8::has_affine_map() const
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



std::unique_ptr<Elem> Hex8::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Quad4,Hex8>>(this,i);

  else
    {
      std::unique_ptr<Elem> face = libmesh_make_unique<Quad4>();
      face->subdomain_id() = this->subdomain_id();

      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(Hex8::side_nodes_map[i][n]);

      return face;
    }
}



std::unique_ptr<Elem> Hex8::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2,Hex8>>(this,i);
}



void Hex8::connectivity(const unsigned int libmesh_dbg_var(sc),
                        const IOPackage iop,
                        std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
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



#ifdef LIBMESH_ENABLE_AMR

const float Hex8::_embedding_matrix[8][8][8] =
  {
    // The 8 children of the Hex-type elements can be thought of as being
    // associated with the 8 vertices of the Hex.  Some of the children are
    // numbered the same as their corresponding vertex, while some are
    // not.  The children which are numbered differently have been marked
    // with ** in the comments below.

    // embedding matrix for child 0 (child 0 is associated with vertex 0)
    {
      //  0     1     2     3     4     5     6     7
      { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
      { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 4
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 5
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 6
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}  // 7
    },

    // embedding matrix for child 1 (child 1 is associated with vertex 1)
    {
      //  0     1     2     3     4     5     6     7
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 3
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 4
      { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 5
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 6
      {.125, .125, .125, .125, .125, .125, .125, .125}  // 7
    },

    // embedding matrix for child 2 (child 2 is associated with vertex 3**)
    {
      //  0      1    2     3     4     5     6     7
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 3
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 4
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 5
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 6
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}  // 7
    },

    // embedding matrix for child 3 (child 3 is associated with vertex 2**)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 4
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 5
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 6
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}  // 7
    },

    // embedding matrix for child 4 (child 4 is associated with vertex 4)
    {
      //  0      1    2     3     4     5     6     7
      { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 1
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 2
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 3
      { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 7
    },

    // embedding matrix for child 5 (child 5 is associated with vertex 5)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 1
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 2
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 6
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}  // 7
    },

    // embedding matrix for child 6 (child 6 is associated with vertex 7**)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 0
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 1
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 7
    },

    // embedding matrix for child 7 (child 7 is associated with vertex 6**)
    {
      //  0      1    2     3     4     5     6     7
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 0
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 1
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 2
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 3
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 7
    }
  };




#endif



Real Hex8::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1), x2 = point(2), x3 = point(3),
    x4 = point(4), x5 = point(5), x6 = point(6), x7 = point(7);

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

} // namespace libMesh
