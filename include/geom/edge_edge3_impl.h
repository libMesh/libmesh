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

#ifndef LIBMESH_EDGE_EDGE3_IMPL_H
#define LIBMESH_EDGE_EDGE3_IMPL_H

// Local includes
#include "libmesh/edge_edge3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

template <typename RealType>
bool Edge3Templ<RealType>::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

template <typename RealType>
bool Edge3Templ<RealType>::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}

template <typename RealType>
bool Edge3Templ<RealType>::is_face(const unsigned int ) const
{
  return false;
}

template <typename RealType>
bool Edge3Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, Edge3::num_nodes);
  return (s == n);
}

template <typename RealType>
bool Edge3Templ<RealType>::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



template <typename RealType>
bool Edge3Templ<RealType>::has_affine_map() const
{
  return (this->point(2).relative_fuzzy_equals
          ((this->point(0) + this->point(1))/2));
}



template <typename RealType>
Order Edge3Templ<RealType>::default_order() const
{
  return SECOND;
}



template <typename RealType>
void Edge3Templ<RealType>::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less_equal (sc, 1);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
        switch (sc)
          {
          case 0:
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(2)+1;
            return;

          case 1:
            conn[0] = this->node_id(2)+1;
            conn[1] = this->node_id(1)+1;
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }


    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;

        /*
          switch (sc)
          {
          case 0:
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(2);

          return;

          case 1:
          conn[0] = this->node_id(2);
          conn[1] = this->node_id(1);

          return;

          default:
          libmesh_error_msg("Invalid sc = " << sc);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



template <typename RealType>
std::pair<unsigned short int, unsigned short int>
Edge3Templ<RealType>::second_order_child_vertex (const unsigned int) const
{
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



template <typename RealType>
RealType Edge3Templ<RealType>::volume () const
{
  // Finding the (exact) length of a general quadratic element
  // is a surprisingly complicated formula.
  Point A = this->point(0) + this->point(1) - 2*this->point(2);
  Point B = (this->point(1) - this->point(0))/2;

  const Real a = A.norm_sq();
  const Real c = B.norm_sq();

  // Degenerate straight line case
  if (a < TOLERANCE*TOLERANCE)
    return 2. * std::sqrt(c);

  const Real b = 2.*(A*B);
  const Real ba=b/a;
  const Real ca=c/a;

  libmesh_assert (1.-ba+ca>0.);

  const Real s1 = std::sqrt(1. - ba + ca);
  const Real s2 = std::sqrt(1. + ba + ca);

  Real log_term = (1. - 0.5*ba + s1) / (-1. - 0.5*ba + s2);
  libmesh_assert(!libmesh_isnan(log_term) && log_term > 0.);

  return 0.5*std::sqrt(a)*((1.-0.5*ba)*s1 +
                           (1.+0.5*ba)*s2 +
                           (ca - 0.25*ba*ba)*std::log(log_term)
                           );
}



template <typename RealType>
BoundingBoxTempl<RealType> Edge3Templ<RealType>::loose_bounding_box () const
{
  // This might be a curved line through 2-space or 3-space, in which
  // case the full bounding box can be larger than the bounding box of
  // just the nodes.
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      Real center = this->point(2)(d);
      Real hd = std::max(std::abs(center - this->point(0)(d)),
                         std::abs(center - this->point(1)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



template <typename RealType>
dof_id_type Edge3Templ<RealType>::key () const
{
  return this->compute_key(this->node_id(2));
}


} // namespace libMesh

#endif // LIBMESH_EDGE_EDGE3_IMPL_H
