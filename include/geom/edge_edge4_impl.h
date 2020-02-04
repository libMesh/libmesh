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

#ifndef LIBMESH_EDGE_EDGE4_IMPL_H
#define LIBMESH_EDGE_EDGE4_IMPL_H

// Local includes
#include "libmesh/edge_edge4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

template <typename RealType>
bool Edge4Templ<RealType>::is_vertex(const unsigned int i) const
{
  return (i==0) || (i==1);
}

template <typename RealType>
bool Edge4Templ<RealType>::is_edge(const unsigned int i) const
{
  return (i==2) || (i==3);
}

template <typename RealType>
bool Edge4Templ<RealType>::is_face(const unsigned int ) const
{
  return false;
}

template <typename RealType>
bool Edge4Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, Edge4::num_nodes);
  return (s == n);
}

template <typename RealType>
bool Edge4Templ<RealType>::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



template <typename RealType>
bool Edge4Templ<RealType>::has_affine_map() const
{
  if (!this->point(2).relative_fuzzy_equals
      ((this->point(0)*2. + this->point(1))/3.))
    return false;
  if (!this->point(3).relative_fuzzy_equals
      ((this->point(0) + this->point(1)*2.)/3.))
    return false;
  return true;
}



template <typename RealType>
Order Edge4Templ<RealType>::default_order() const
{
  return THIRD;
}



template <typename RealType>
void Edge4Templ<RealType>::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less_equal (sc, 2);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch(iop)
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
            conn[1] = this->node_id(3)+1;
            return;

          case 2:
            conn[0] = this->node_id(3)+1;
            conn[1] = this->node_id(1)+1;
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }

      }

    case VTK:
      {

        switch (sc)
          {
          case 0:
            conn[0] = this->node_id(0);
            conn[1] = this->node_id(2);
            return;

          case 1:
            conn[0] = this->node_id(2);
            conn[1] = this->node_id(3);
            return;

          case 2:
            conn[0] = this->node_id(3);
            conn[1] = this->node_id(1);
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



template <typename RealType>
BoundingBoxTempl<RealType> Edge4Templ<RealType>::loose_bounding_box () const
{
  // This might be a curved line through 2-space or 3-space, in which
  // case the full bounding box can be larger than the bounding box of
  // just the nodes.
  //
  // FIXME - I haven't yet proven the formula below to be correct for
  // cubics - RHS
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      Real center = (this->point(2)(d) + this->point(3)(d))/2;
      Real hd = std::max(std::abs(center - this->point(0)(d)),
                         std::abs(center - this->point(1)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



template <typename RealType>
dof_id_type Edge4Templ<RealType>::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2),
                           this->node_id(3));
}



template <typename RealType>
RealType Edge4Templ<RealType>::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0),
    x1 = this->point(1),
    x2 = this->point(2),
    x3 = this->point(3);

  // Construct constant data vectors.
  // \vec{x}_{\xi}  = \vec{a1}*xi**2 + \vec{b1}*xi + \vec{c1}
  // This is copy-pasted directly from the output of a Python script.
  Point
    a1 = -27*x0/16 + 27*x1/16 + 81*x2/16 - 81*x3/16,
    b1 = 9*x0/8 + 9*x1/8 - 9*x2/8 - 9*x3/8,
    c1 = x0/16 - x1/16 - 27*x2/16 + 27*x3/16;

  // 4 point quadrature, 7th-order accurate
  const unsigned int N = 4;
  const Real q[N] = {-std::sqrt(525 + 70*std::sqrt(30.)) / 35,
                     -std::sqrt(525 - 70*std::sqrt(30.)) / 35,
                     std::sqrt(525 - 70*std::sqrt(30.)) / 35,
                     std::sqrt(525 + 70*std::sqrt(30.)) / 35};
  const Real w[N] = {(18 - std::sqrt(30.)) / 36,
                     (18 + std::sqrt(30.)) / 36,
                     (18 + std::sqrt(30.)) / 36,
                     (18 - std::sqrt(30.)) / 36};

  Real vol=0.;
  for (unsigned int i=0; i<N; ++i)
    vol += w[i] * (q[i]*q[i]*a1 + q[i]*b1 + c1).norm();

  return vol;
}

} // namespace libMesh

#endif // LIBMESH_EDGE_EDGE4_IMPL_H
