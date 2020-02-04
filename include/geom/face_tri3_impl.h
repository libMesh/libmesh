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

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Tri3 class member functions

template <typename RealType>
bool Tri3Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Tri3Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Tri3Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Tri3Templ<RealType>::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Tri3Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

template <typename RealType>
Order Tri3Templ<RealType>::default_order() const
{
  return FIRST;
}

template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Tri3Templ<RealType>::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Edge2,Tri3>>(this,i);

  else
    {
      std::unique_ptr<Elem> edge = libmesh_make_unique<Edge2>();
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : edge->node_index_range())
        edge->set_node(n) = this->node_ptr(Tri3::side_nodes_map[i][n]);

      return edge;
    }
}



template <typename RealType>
void Tri3Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->template simple_build_side_ptr<Tri3>(side, i, EDGE2);
}



template <typename RealType>
void Tri3Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sf),
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
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        return;
      }

    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}






template <typename RealType>
RealType Tri3Templ<RealType>::volume () const
{
  // 3-node triangles have the following formula for computing the area
  return 0.5 * cross_norm(this->point(1) - this->point(0),
                          this->point(2) - this->point(0));
}



template <typename RealType>
std::pair<Real, Real> Tri3Templ<RealType>::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
    len_10=v10.norm(),
    len_20=v20.norm(),
    len_21=v21.norm();

  const Real
    theta0=std::acos(( v10*v20)/len_10/len_20),
    theta1=std::acos((-v10*v21)/len_10/len_21),
    theta2=libMesh::pi - theta0 - theta1
    ;

  libmesh_assert_greater (theta0, 0.);
  libmesh_assert_greater (theta1, 0.);
  libmesh_assert_greater (theta2, 0.);

  return std::make_pair(std::min(theta0, std::min(theta1,theta2)),
                        std::max(theta0, std::max(theta1,theta2)));
}

template <typename RealType>
bool Tri3Templ<RealType>::contains_point (const Point & p, Real tol) const
{
  // See "Barycentric Technique" section at
  // http://www.blackpawn.com/texts/pointinpoly for details.

  // Compute vectors
  Point v0 = this->point(1) - this->point(0);
  Point v1 = this->point(2) - this->point(0);
  Point v2 = p - this->point(0);

  // Compute dot products
  Real dot00 = v0 * v0;
  Real dot01 = v0 * v1;
  Real dot02 = v0 * v2;
  Real dot11 = v1 * v1;
  Real dot12 = v1 * v2;

  // Out of plane check
  if (std::abs(triple_product(v2, v0, v1)) / std::max(dot00, dot11) > tol)
    return false;

  // Compute barycentric coordinates
  Real invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  Real u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  Real v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  return (u > -tol) && (v > -tol) && (u + v < 1 + tol);
}

template <typename RealType>
BoundingBoxTempl<RealType>
Tri3Templ<RealType>::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


} // namespace libMesh
