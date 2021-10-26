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
#include "libmesh/cell_hex8.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/fe_lagrange_shape_1D.h"

// C++ includes
#include <array>

namespace libMesh
{




// ------------------------------------------------------------
// Hex8 class static member initializations
const int Hex8::num_nodes;
const int Hex8::num_sides;
const int Hex8::num_edges;
const int Hex8::num_children;
const int Hex8::nodes_per_side;
const int Hex8::nodes_per_edge;

const unsigned int Hex8::side_nodes_map[Hex8::num_sides][Hex8::nodes_per_side] =
  {
    {0, 3, 2, 1}, // Side 0
    {0, 1, 5, 4}, // Side 1
    {1, 2, 6, 5}, // Side 2
    {2, 3, 7, 6}, // Side 3
    {3, 0, 4, 7}, // Side 4
    {4, 5, 6, 7}  // Side 5
  };

const unsigned int Hex8::edge_nodes_map[Hex8::num_edges][Hex8::nodes_per_edge] =
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

const unsigned int Hex8::edge_sides_map[Hex8::num_edges][2] =
  {
    {0, 1}, // Edge 0
    {0, 2}, // Edge 1
    {0, 3}, // Edge 2
    {0, 4}, // Edge 3
    {1, 4}, // Edge 4
    {1, 2}, // Edge 5
    {2, 3}, // Edge 6
    {3, 4}, // Edge 7
    {1, 5}, // Edge 8
    {2, 5}, // Edge 9
    {3, 5}, // Edge 10
    {4, 5}  // Edge 11
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
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Hex8::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Hex8::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Hex8::is_node_on_edge(const unsigned int n,
                           const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



bool Hex8::has_affine_map() const
{
  // Make sure x-edge endpoints are affine
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(4), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(7), affine_tol))
    return false;
  // Make sure xz-faces are identical parallelograms
  v = this->point(4) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(3), affine_tol))
    return false;
  // If all the above checks out, the map is affine
  return true;
}



Order Hex8::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Hex8::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  return this->simple_build_side_ptr<Quad4, Hex8>(i, proxy);
}



void Hex8::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->simple_build_side_ptr<Hex8>(side, i, QUAD4);
}



std::unique_ptr<Elem> Hex8::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge2,Hex8>(i);
}



void Hex8::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Hex8>(edge, i, EDGE2);
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

const Real Hex8::_embedding_matrix[Hex8::num_children][Hex8::num_nodes][Hex8::num_nodes] =
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



Point Hex8::centroid_from_points(
  const Point & x0, const Point & x1, const Point & x2, const Point & x3,
  const Point & x4, const Point & x5, const Point & x6, const Point & x7)
{
  // The Jacobian is dx/d(xi) dot (dx/d(eta) cross dx/d(zeta)), where
  // dx/d(xi)   = a1*eta*zeta + b1*eta + c1*zeta + d1
  // dx/d(eta)  = a2*xi*zeta  + b2*xi  + c2*zeta + d2
  // dx/d(zeta) = a3*xi*eta   + b3*xi  + c3*eta  + d3

  // Notes:
  // 1.) Several of these coefficient vectors are equal, as noted below.
  // 2.) These are all off by a factor of 8, but this cancels when we
  //     divide by the volume, which will also be off by the same
  //     factor.
  Point
    a1 = -x0 + x1 - x2 + x3 + x4 - x5 + x6 - x7,
    a2 = a1,
    a3 = a1;

  Point
    b1 = x0 - x1 + x2 - x3 + x4 - x5 + x6 - x7,
    b2 = b1,
    b3 = x0 - x1 - x2 + x3 - x4 + x5 + x6 - x7;

  Point
    c1 = b3,
    c2 = x0 + x1 - x2 - x3 - x4 - x5 + x6 + x7,
    c3 = c2;

  Point
    d1 = -x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7,
    d2 = -x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7,
    d3 = -x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7;

  // Use 2x2x2 quadrature to compute the integral of each basis
  // function (as defined on the [-1,1]^3 reference domain). We use
  // a quadrature rule which is exact for tri-cubics. The weights for
  // this rule are all equal to 1.
  static const Real q[2] = {-std::sqrt(3.)/3, std::sqrt(3.)/3.};

  // Indices for computing tensor product basis functions. This is
  // copied from fe_lagrange_shape_3D.C
  static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
  static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
  static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

  // Compute nodal volumes
  std::array<Real, Hex8::num_nodes> V{};

  for (const auto & xi : q)
    for (const auto & eta : q)
      for (const auto & zeta : q)
      {
        Real jxw = triple_product(a1*eta*zeta + b1*eta + c1*zeta + d1,
                                  a2*xi*zeta  + b2*xi  + c2*zeta + d2,
                                  a3*xi*eta   + b3*xi  + c3*eta  + d3);

        for (int i=0; i<Hex8::num_nodes; ++i)
          V[i] += jxw *
            fe_lagrange_1D_linear_shape(i0[i], xi) *
            fe_lagrange_1D_linear_shape(i1[i], eta) *
            fe_lagrange_1D_linear_shape(i2[i], zeta);
      }

  // Compute centroid
  return
    (x0*V[0] + x1*V[1] + x2*V[2] + x3*V[3] + x4*V[4] + x5*V[5] + x6*V[6] + x7*V[7]) /
    (V[0] + V[1] + V[2] + V[3] + V[4] + V[5] + V[6] + V[7]);
}


Point Hex8::true_centroid () const
{
  return Hex8::centroid_from_points
    (point(0), point(1), point(2), point(3),
     point(4), point(5), point(6), point(7));
}


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

BoundingBox
Hex8::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


void
Hex8::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 24);
  const unsigned int side = perm_num % 6;
  const unsigned int rotate = perm_num / 6;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(4,5,6,7);
      swap4neighbors(1,2,3,4);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap4nodes(3,7,4,0);
    swap4nodes(2,6,5,1);
    swap4neighbors(0,3,5,1);
    break;
  case 2:
    swap4nodes(0,4,5,1);
    swap4nodes(3,7,6,2);
    swap4neighbors(0,4,5,2);
    break;
  case 3:
    swap4nodes(0,4,7,3);
    swap4nodes(1,5,6,2);
    swap4neighbors(0,1,5,3);
    break;
  case 4:
    swap4nodes(1,5,4,0);
    swap4nodes(2,6,7,3);
    swap4neighbors(0,2,5,4);
    break;
  case 5:
    swap2nodes(0,7);
    swap2nodes(1,6);
    swap2nodes(2,5);
    swap2nodes(3,4);
    swap2neighbors(0,5);
    swap2neighbors(1,3);
    break;
  default:
    libmesh_error();
  }
}


} // namespace libMesh
