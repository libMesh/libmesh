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

#ifndef LIBMESH_FACE_TRI_IMPL_H
#define LIBMESH_FACE_TRI_IMPL_H

// Local includes
#include "libmesh/face_tri.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_elem_quality.h"

namespace libMesh
{

// ------------------------------------------------------------
// Tri class member functions
template <typename RealType>
dof_id_type TriTempl<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tri3::side_nodes_map[s][0]),
                           this->node_id(Tri3::side_nodes_map[s][1]));
}



template <typename RealType>
unsigned int TriTempl<RealType>::which_node_am_i(unsigned int side,
                                  unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 2);

  return Tri3::side_nodes_map[side][side_node];
}



template <typename RealType>
dof_id_type TriTempl<RealType>::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2));
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> TriTempl<RealType>::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> edge = libmesh_make_unique<Edge2>();

  for (auto n : edge->node_index_range())
    edge->set_node(n) = this->node_ptr(Tri3::side_nodes_map[i][n]);

  return edge;
}



template <typename RealType>
void TriTempl<RealType>::side_ptr (std::unique_ptr<Elem> & side,
                    const unsigned int i)
{
  this->template simple_side_ptr<Tri,Tri3>(side, i, EDGE2);
}



template <typename RealType>
bool TriTempl<RealType>::is_child_on_side(const unsigned int c,
                           const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (c == s || c == (s+1)%3);
}



template <typename RealType>
Real TriTempl<RealType>::quality (const ElemQuality q) const
{
  switch (q)
    {

      /**
       * Source: Netgen, meshtool.cpp, TriangleQualityInst
       */
    case DISTORTION:
    case STRETCH:
      {
        const Point & p1 = this->point(0);
        const Point & p2 = this->point(1);
        const Point & p3 = this->point(2);

        Point v1 = p2 - p1;
        Point v2 = p3 - p1;
        Point v3 = p3 - p2;
        const Real l1 = v1.norm();
        const Real l2 = v2.norm();
        const Real l3 = v3.norm();

        // if one length is 0, quality is quite bad!
        if ((l1 <=0.) || (l2 <= 0.) || (l3 <= 0.))
          return 0.;

        const Real s1 = std::sin(std::acos(v1*v2/l1/l2)/2.);
        v1 *= -1;
        const Real s2 = std::sin(std::acos(v1*v3/l1/l3)/2.);
        const Real s3 = std::sin(std::acos(v2*v3/l2/l3)/2.);

        return 8. * s1 * s2 * s3;

      }
    default:
      return Elem::quality(q);
    }

  /**
   * I don't know what to do for this metric.
   * Maybe the base class knows.  We won't get
   * here because of the default case above.
   */
  return Elem::quality(q);

}






template <typename RealType>
std::pair<Real, Real> TriTempl<RealType>::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case MAX_ANGLE:
      bounds.first  = 60.;
      bounds.second = 90.;
      break;

    case MIN_ANGLE:
      bounds.first  = 30.;
      bounds.second = 60.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 1.3;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.155;
      break;

    case SIZE:
    case SHAPE:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}

} // namespace libMesh

#endif // LIBMESH_FACE_TRI_IMPL_H
