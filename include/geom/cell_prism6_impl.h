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

#ifndef LIBMESH_CELL_PRISM6_IMPL_H
#define LIBMESH_CELL_PRISM6_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Prism6 class member functions

template <typename RealType>
bool Prism6Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Prism6Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Prism6Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Prism6Templ<RealType>::is_node_on_side(const unsigned int n,
                             const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Prism6Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  auto trim = (s > 0 && s < 4) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

template <typename RealType>
bool Prism6Templ<RealType>::is_node_on_edge(const unsigned int n,
                             const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



template <typename RealType>
bool Prism6Templ<RealType>::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2)))
    return false;
  return true;
}



template <typename RealType>
Order Prism6Templ<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Prism6Templ<RealType>::build_side_ptr (const unsigned int i,
                                              bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch(i)
        {
        case 0:
        case 4:
          return libmesh_make_unique<Side<Tri3,Prism6>>(this,i);

        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<Quad4,Prism6>>(this,i);

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
        case 0: // the triangular face at z=-1
        case 4: // the triangular face at z=1
          {
            face = libmesh_make_unique<Tri3>();
            break;
          }
        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = libmesh_make_unique<Quad4>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Prism6::side_nodes_map[i][n]);

      return face;
    }
}



template <typename RealType>
void Prism6Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                             const unsigned int i)
{
  this->side_ptr(side, i);
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Prism6Templ<RealType>::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2,Prism6>>(this,i);
}



template <typename RealType>
void Prism6Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[3] = this->node_id(2)+1;
        conn[4] = this->node_id(3)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(5)+1;
        conn[7] = this->node_id(5)+1;
        return;
      }

    case VTK:
      {
        conn.resize(6);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(2);
        conn[2] = this->node_id(1);
        conn[3] = this->node_id(3);
        conn[4] = this->node_id(5);
        conn[5] = this->node_id(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


template <typename RealType>
RealType Prism6Templ<RealType>::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0), x1 = this->point(1), x2 = this->point(2),
    x3 = this->point(3), x4 = this->point(4), x5 = this->point(5);

  // constant and zeta terms only.  These are copied directly from a
  // Python script.
  Point dx_dxi[2] =
    {
      -x0/2 + x1/2 - x3/2 + x4/2, // constant
      x0/2 - x1/2 - x3/2 + x4/2,  // zeta
    };

  // constant and zeta terms only.  These are copied directly from a
  // Python script.
  Point dx_deta[2] =
    {
      -x0/2 + x2/2 - x3/2 + x5/2, // constant
      x0/2 - x2/2 - x3/2 + x5/2,  // zeta
    };

  // Constant, xi, and eta terms
  Point dx_dzeta[3] =
    {
      -x0/2 + x3/2,              // constant
      x0/2 - x2/2 - x3/2 + x5/2, // eta
      x0/2 - x1/2 - x3/2 + x4/2  // xi
    };

  // The quadrature rule the Prism6 is a tensor product between a
  // four-point TRI3 rule (in xi, eta) and a two-point EDGE2 rule (in
  // zeta) which is capable of integrating cubics exactly.

  // Number of points in the 2D quadrature rule.
  const int N2D = 4;

  static const Real w2D[N2D] =
    {
      Real(1.5902069087198858469718450103758e-01L),
      Real(9.0979309128011415302815498962418e-02L),
      Real(1.5902069087198858469718450103758e-01L),
      Real(9.0979309128011415302815498962418e-02L)
    };

  static const Real xi[N2D] =
    {
      Real(1.5505102572168219018027159252941e-01L),
      Real(6.4494897427831780981972840747059e-01L),
      Real(1.5505102572168219018027159252941e-01L),
      Real(6.4494897427831780981972840747059e-01L)
    };

  static const Real eta[N2D] =
    {
      Real(1.7855872826361642311703513337422e-01L),
      Real(7.5031110222608118177475598324603e-02L),
      Real(6.6639024601470138670269327409637e-01L),
      Real(2.8001991549907407200279599420481e-01L)
    };

  // Number of points in the 1D quadrature rule.  The weights of the
  // 1D quadrature rule are equal to 1.
  const int N1D = 2;

  // Points of the 1D quadrature rule
  static const Real zeta[N1D] =
    {
      -std::sqrt(3.)/3,
      std::sqrt(3.)/3.
    };

  RealType vol = 0.;
  for (int i=0; i<N2D; ++i)
    {
      // dx_dzeta depends only on the 2D quadrature rule points.
      Point dx_dzeta_q = dx_dzeta[0] + eta[i]*dx_dzeta[1] + xi[i]*dx_dzeta[2];

      for (int j=0; j<N1D; ++j)
        {
          // dx_dxi and dx_deta only depend on the 1D quadrature rule points.
          Point
            dx_dxi_q  = dx_dxi[0]  + zeta[j]*dx_dxi[1],
            dx_deta_q = dx_deta[0] + zeta[j]*dx_deta[1];

          // Compute scalar triple product, multiply by weight, and accumulate volume.
          vol += w2D[i] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q);
        }
    }

  return vol;
}

template <typename RealType>
BoundingBoxTempl<RealType>
Prism6Templ<RealType>::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM6_IMPL_H
