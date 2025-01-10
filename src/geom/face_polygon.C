// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/face_polygon.h"

// Local includes
#include "libmesh/edge_edge2.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/hashword.h"

// C++ includes
#include <array>


namespace libMesh
{

// ------------------------------------------------------------
// Polygon class static member initializations
const int Polygon::num_children;

// ------------------------------------------------------------
// Polygon class member functions
Point Polygon::master_point (const unsigned int i) const
{
  const unsigned int ns = this->n_sides();

  // Center-to-midside to center-to-vertex angle
  const Real pi_over_ns = libMesh::pi / ns;

  // Center-to-midside distance
  const Real min_r = 0.5/tan(pi_over_ns);

  // Center-to-vertex distance
  const Real max_r = std::sqrt(min_r*min_r + 0.25);

  const Point center(0.5, min_r);

  libmesh_assert_less(i, this->n_nodes());

  return center + Point(max_r*sin((int(i)*2-1)*pi_over_ns),
                        -max_r*cos((int(i)*2-1)*pi_over_ns));

}


dof_id_type Polygon::key (const unsigned int s) const
{
  const int ns = this->n_sides();
  libmesh_assert_less (s, ns);

  return this->compute_key(this->node_id(s),
                           this->node_id((s+1)%ns));
}



dof_id_type Polygon::low_order_key (const unsigned int s) const
{
  const int ns = this->n_sides();
  libmesh_assert_less (s, ns);

  return this->compute_key(this->node_id(s),
                           this->node_id((s+1)%ns));
}



unsigned int Polygon::local_side_node(unsigned int side,
                                      unsigned int side_node) const
{
  const auto ns = this->n_sides();
  libmesh_assert_less (side, ns);
  libmesh_assert_less (side_node, this->n_nodes()/ns + 1);

  if (!side_node)
    return side;

  if (side_node == 1)
    return (side + 1) % ns;

  return (side_node - 1) * ns + side;
}



unsigned int Polygon::local_edge_node(unsigned int edge,
                                      unsigned int edge_node) const
{
  return local_side_node(edge, edge_node);
}



dof_id_type Polygon::key () const
{
  std::vector<dof_id_type> node_ids;
  for (const auto & n : this->node_ref_range())
    node_ids.push_back(n.id());

  return Utility::hashword(node_ids);
}



std::unique_ptr<Elem> Polygon::side_ptr (const unsigned int i)
{
  const auto ns = this->n_sides();
  libmesh_assert_less (i, ns);

  std::unique_ptr<Elem> edge = std::make_unique<Edge2>();

  edge->set_node(0) = this->node_ptr(i);
  edge->set_node(1) = this->node_ptr((i+1)%ns);

  return edge;
}



void Polygon::side_ptr (std::unique_ptr<Elem> & side,
                        const unsigned int i)
{
  const auto ns = this->n_sides();
  libmesh_assert_less (i, ns);

  if (!side.get() || side->type() != EDGE2)
    side = std::make_unique<Edge2>();

  side->subdomain_id() = this->subdomain_id();

  side->set_node(0) = this->node_ptr(i);
  side->set_node(1) = this->node_ptr((i+1)%ns);
}



bool Polygon::is_child_on_side(const unsigned int /*c*/,
                               const unsigned int /*s*/) const
{
  libmesh_not_implemented();
  return false;
}



unsigned int Polygon::opposite_side(const unsigned int side_in) const
{
  const auto ns = this->n_sides();
  if (ns % 2)
    libmesh_error();

  return (ns / 2 + side_in) % ns;
}



bool Polygon::is_flipped() const
{
  // Just check based on vertices, the first ns nodes
  const auto ns = this->n_sides();

#if LIBMESH_DIM > 2
  // Don't bother outside the XY plane
  for (auto n : make_range(ns))
    if (this->point(n)(2))
      return false;
#endif

  Real twice_signed_area = 0;
  for (auto n : make_range(ns))
    {
      const Point & p1 = this->point(n);
      const Point & p2 = this->point((n+1)%ns);
      twice_signed_area += p1(0)*p2(1)-p1(1)*p2(0);
    }

  return (twice_signed_area < 0);
}


std::vector<unsigned int>
Polygon::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  const unsigned int ns = this->n_sides();

  // For vertices, we have two adjacent edges; otherwise each of the
  // mid-edge nodes is adjacent only to the edge it is on.
  if (this->is_vertex(n))
    return {n, (n+ns-1)%ns};

  return {n % ns};
}


std::pair<Real, Real> Polygon::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {
    case EDGE_LENGTH_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case MIN_ANGLE:
      bounds.first  = 90. - (180./this->n_sides());
      bounds.second = 180. - (360./this->n_sides());
      break;

    case MAX_ANGLE:
      bounds.first  = 180. - (360./this->n_sides());
      bounds.second = 180. - (180./this->n_sides());
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}


std::array<Point, 3> Polygon::master_subtriangle (unsigned int i) const
{
  libmesh_assert_less(i, this->_triangulation.size());

  const auto & tri = this->_triangulation[i];

  return { this->master_point(tri[0]),
           this->master_point(tri[1]),
           this->master_point(tri[2]) };
}


} // namespace libMesh
