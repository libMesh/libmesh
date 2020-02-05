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

#ifndef LIBMESH_FACE_QUAD9_IMPL_H
#define LIBMESH_FACE_QUAD9_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Quad9 class member functions

template <typename RealType>
bool Quad9Templ<RealType>::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

template <typename RealType>
bool Quad9Templ<RealType>::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  if (i > 7)
    return false;
  return true;
}

template <typename RealType>
bool Quad9Templ<RealType>::is_face(const unsigned int i) const
{
  if (i > 7)
    return true;
  return false;
}

template <typename RealType>
bool Quad9Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Quad9Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

template <typename RealType>
bool Quad9Templ<RealType>::has_affine_map() const
{
  // make sure corners form a parallelogram
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3)))
    return false;
  // make sure "horizontal" sides are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(3)))
    return false;
  // make sure "vertical" sides are straight
  // and the center node is centered
  v = (this->point(3) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(8) - this->point(4)))
    return false;
  return true;
}



template <typename RealType>
Order Quad9Templ<RealType>::default_order() const
{
  return SECOND;
}



template <typename RealType>
dof_id_type Quad9Templ<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:

      return
        this->compute_key (this->node_id(4));

    case 1:

      return
        this->compute_key (this->node_id(5));

    case 2:

      return
        this->compute_key (this->node_id(6));

    case 3:

      return
        this->compute_key (this->node_id(7));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



template <typename RealType>
dof_id_type Quad9Templ<RealType>::key () const
{
  return this->compute_key(this->node_id(8));
}



template <typename RealType>
unsigned int Quad9Templ<RealType>::which_node_am_i(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Quad9::nodes_per_side);

  return Quad9::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Quad9Templ<RealType>::build_side_ptr (const unsigned int i,
                                             bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Edge3,Quad9>>(this,i);

  else
    {
      std::unique_ptr<Elem> edge = libmesh_make_unique<Edge3>();
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : edge->node_index_range())
        edge->set_node(n) = this->node_ptr(Quad9::side_nodes_map[i][n]);

      return edge;
    }
}



template <typename RealType>
void Quad9Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->template simple_build_side_ptr<Quad9>(side, i, EDGE3);
}



template <typename RealType>
void Quad9Templ<RealType>::connectivity(const unsigned int sf,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        switch(sf)
          {
          case 0:
            // linear sub-quad 0
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(8)+1;
            conn[3] = this->node_id(7)+1;
            return;

          case 1:
            // linear sub-quad 1
            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(8)+1;
            return;

          case 2:
            // linear sub-quad 2
            conn[0] = this->node_id(7)+1;
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(3)+1;
            return;

          case 3:
            // linear sub-quad 3
            conn[0] = this->node_id(8)+1;
            conn[1] = this->node_id(5)+1;
            conn[2] = this->node_id(2)+1;
            conn[3] = this->node_id(6)+1;
            return;

          default:
            libmesh_error_msg("Invalid sf = " << sf);
          }
      }

    case VTK:
      {
        conn.resize(9);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        conn[3] = this->node_id(3);
        conn[4] = this->node_id(4);
        conn[5] = this->node_id(5);
        conn[6] = this->node_id(6);
        conn[7] = this->node_id(7);
        conn[8] = this->node_id(8);
        return;

        /*
          switch(sf)
          {
          case 0:
          // linear sub-quad 0
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(4);
          conn[2] = this->node_id(8);
          conn[3] = this->node_id(7);

          return;

          case 1:
          // linear sub-quad 1
          conn[0] = this->node_id(4);
          conn[1] = this->node_id(1);
          conn[2] = this->node_id(5);
          conn[3] = this->node_id(8);

          return;

          case 2:
          // linear sub-quad 2
          conn[0] = this->node_id(7);
          conn[1] = this->node_id(8);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(3);

          return;

          case 3:
          // linear sub-quad 3
          conn[0] = this->node_id(8);
          conn[1] = this->node_id(5);
          conn[2] = this->node_id(2);
          conn[3] = this->node_id(6);

          return;

          default:
          libmesh_error_msg("Invalid sf = " << sf);
          }*/
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



template <typename RealType>
BoundingBoxTempl<RealType> Quad9Templ<RealType>::loose_bounding_box () const
{
  // This might have curved edges, or might be a curved surface in
  // 3-space, in which case the full bounding box can be larger than
  // the bounding box of just the nodes.
  //
  //
  // FIXME - I haven't yet proven the formula below to be correct for
  // biquadratics - RHS
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      const Real center = this->point(8)(d);
      Real hd = std::abs(center - this->point(0)(d));
      for (unsigned int p=0; p != 8; ++p)
        hd = std::max(hd, std::abs(center - this->point(p)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



template <typename RealType>
RealType Quad9Templ<RealType>::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0), x1 = this->point(1), x2 = this->point(2),
    x3 = this->point(3), x4 = this->point(4), x5 = this->point(5),
    x6 = this->point(6), x7 = this->point(7), x8 = this->point(8);

  // Construct constant data vectors.
  // \vec{x}_{\xi}  = \vec{a1}*xi*eta^2 + \vec{b1}*eta**2 + \vec{c1}*xi*eta + \vec{d1}*xi + \vec{e1}*eta + \vec{f1}
  // \vec{x}_{\eta} = \vec{a2}*xi^2*eta + \vec{b2}*xi**2  + \vec{c2}*xi*eta + \vec{d2}*xi + \vec{e2}*eta + \vec{f2}
  // This is copy-pasted directly from the output of a Python script.
  Point
    a1 = x0/2 + x1/2 + x2/2 + x3/2 - x4 - x5 - x6 - x7 + 2*x8,
    b1 = -x0/4 + x1/4 + x2/4 - x3/4 - x5/2 + x7/2,
    c1 = -x0/2 - x1/2 + x2/2 + x3/2 + x4 - x6,
    d1 = x5 + x7 - 2*x8,
    e1 = x0/4 - x1/4 + x2/4 - x3/4,
    f1 = x5/2 - x7/2,
    a2 = a1,
    b2 = -x0/4 - x1/4 + x2/4 + x3/4 + x4/2 - x6/2,
    c2 = -x0/2 + x1/2 + x2/2 - x3/2 - x5 + x7,
    d2 = x0/4 - x1/4 + x2/4 - x3/4,
    e2 = x4 + x6 - 2*x8,
    f2 = -x4/2 + x6/2;

  // 3x3 quadrature, exact for bi-quintics
  const unsigned int N = 3;
  const Real q[N] = {-std::sqrt(15)/5., 0., std::sqrt(15)/5.};
  const Real w[N] = {5./9, 8./9, 5./9};

  RealType vol=0.;
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      vol += w[i] * w[j] *
        cross_norm(q[i]*q[j]*q[j]*a1 + q[j]*q[j]*b1 + q[j]*q[i]*c1 + q[i]*d1 + q[j]*e1 + f1,
                   q[i]*q[i]*q[j]*a2 + q[i]*q[i]*b2 + q[j]*q[i]*c2 + q[i]*d2 + q[j]*e2 + f2);

  return vol;
}




template <typename RealType>
unsigned int Quad9Templ<RealType>::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 4:
    case 5:
    case 6:
    case 7:
      return 2;

    case 8:
      return 4;

    default:
      libmesh_error_msg("Invalid n = " << n);
    }
}



template <typename RealType>
unsigned short int Quad9Templ<RealType>::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
    case 8:
      {
        libmesh_assert_less (v, 4);
        return static_cast<unsigned short int>(v);
      }

    default:
      {
        libmesh_assert_less (v, 2);
        // use the matrix that we inherited from \p Quad
        return this->_second_order_adjacent_vertices[n-this->n_vertices()][v];
      }
    }
}



template <typename RealType>
std::pair<unsigned short int, unsigned short int>
Quad9Templ<RealType>::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in face_quad.C, since they are identical
   * for Quad8 and Quad9 (for the first 4 higher-order nodes)
   */
  return std::pair<unsigned short int, unsigned short int>
    (this->_second_order_vertex_child_number[n],
     this->_second_order_vertex_child_index[n]);
}

} // namespace libMesh

#endif // LIBMESH_FACE_QUAD9_IMPL_H
