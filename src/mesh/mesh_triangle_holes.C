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


#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/mesh_triangle_holes.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/simple_range.h"

namespace
{
  using namespace libMesh;

  // Returns a positive distance iff the ray from source in the
  // direction of ray_target intersects the given edge, -1 otherwise
  Real find_intersection(const Point & source,
                         const Point & ray_target,
                         const Point & edge_pt0,
                         const Point & edge_pt1)
  {
    // Calculate intersection parameters (fractions of the distance
    // along each segment)
    const Real raydx = ray_target(0)-source(0),
               raydy = ray_target(1)-source(1),
               edgedx = edge_pt1(0)-edge_pt0(0),
               edgedy = edge_pt1(1)-edge_pt0(1);
    const Real denom = edgedx * raydy - edgedy * raydx;

    // divide-by-zero means the segments are parallel
    if (denom == 0)
      return -1;

    const Real one_over_denom = 1 / denom;

    const Real targetsdx = edge_pt1(0)-ray_target(0),
               targetsdy = edge_pt1(1)-ray_target(1);

    const Real t_num = targetsdx * raydy -
                       targetsdy * raydx;
    const Real t = t_num * one_over_denom;

    // There's an intersection between the ray line and the edge?
    if (t >= 0 && t < 1)
      {
        const Real u_num = targetsdx * edgedy - targetsdy * edgedx;
        const Real u = u_num * one_over_denom;
        const Real ray_fraction = (1-u);

        // Intersection is in the other direction!?
        if (ray_fraction < 0)
          return -1;

        const Real distance =
          ray_fraction * std::sqrt(raydx*raydx + raydy*raydy);
        return distance;
      }

    return -1;
  }

}


namespace libMesh
{

//
// Hole member functions
//
Real TriangulatorInterface::Hole::area() const
{
  const unsigned int np = this->n_points();

  if (np < 3)
    return 0;

  const Point p0 = this->point(0);

  // Every segment (p_{i-1},p_i) from i=2 on defines a triangle w.r.t.
  // p_0.  Add up the cross products of those triangles.  We'll save
  // the division by 2 and the norm for the end.
  //
  // Your hole points had best be coplanar, but this should work
  // regardless of which plane they're in.

  Point areavec = 0;

  for (unsigned int i=2; i != np; ++i)
    {
      const Point e_0im = this->point(i-1) - p0,
                  e_0i  = this->point(i) - p0;

      areavec += e_0im.cross(e_0i);
    }

  return areavec.norm() / 2;
}


//
// PolygonHole member functions
//
TriangulatorInterface::PolygonHole::PolygonHole(const Point & center,
                                                Real radius,
                                                unsigned int n_points_in) :
  _center(center),
  _radius(radius),
  _n_points(n_points_in)
{}


unsigned int TriangulatorInterface::PolygonHole::n_points() const
{
  return _n_points;
}


Point TriangulatorInterface::PolygonHole::point(const unsigned int n) const
{
  // The nth point lies at the angle theta = 2 * pi * n / _n_points
  const Real theta = static_cast<Real>(n) * 2.0 * libMesh::pi / static_cast<Real>(_n_points);

  return Point(_center(0) + _radius*std::cos(theta), // x=r*cos(theta)
               _center(1) + _radius*std::sin(theta), // y=r*sin(theta)
               0.);
}


Point TriangulatorInterface::PolygonHole::inside() const
{
  // The center of the hole is definitely inside.
  return _center;
}


//
// AffineHole member functions
//
Point TriangulatorInterface::AffineHole::point(const unsigned int n) const
{
  return this->transform(_underlying.point(n));
}


Point TriangulatorInterface::AffineHole::inside() const
{
  return this->transform(_underlying.inside());
}


Point TriangulatorInterface::AffineHole::transform(const Point & p) const
{
  const Real cos_a = std::cos(_angle);
  const Real sin_a = std::sin(_angle);
  return Point(p(0)*cos_a-p(1)*sin_a + _shift(0),
               p(1)*cos_a+p(1)*sin_a + _shift(1));
}


//
// ArbitraryHole member functions
//
TriangulatorInterface::ArbitraryHole::ArbitraryHole(const Point & center,
                                                    std::vector<Point> points)
  : _center(center),
    _points(std::move(points))
{
  _segment_indices.push_back(0);
  _segment_indices.push_back(_points.size());
}


TriangulatorInterface::ArbitraryHole::ArbitraryHole(const Point & center,
                                                    std::vector<Point> points,
                                                    std::vector<unsigned int> segment_indices)
  : _center(center),
    _points(std::move(points)),
    _segment_indices(std::move(segment_indices))
{}


TriangulatorInterface::ArbitraryHole::ArbitraryHole(const Hole & orig)
  : _center(orig.inside())
{
  const unsigned int np = orig.n_points();
  _points.reserve(np);
  for (auto i : make_range(np))
    _points.push_back(orig.point(i));
}


unsigned int TriangulatorInterface::ArbitraryHole::n_points() const
{
  return _points.size();
}


Point TriangulatorInterface::ArbitraryHole::point(const unsigned int n) const
{
  libmesh_assert_less (n, _points.size());
  return _points[n];
}


Point TriangulatorInterface::ArbitraryHole::inside() const
{
  return _center;
}


std::vector<unsigned int> TriangulatorInterface::ArbitraryHole::segment_indices() const
{
  return _segment_indices;
}


//
// MeshedHole member functions
//


TriangulatorInterface::MeshedHole::MeshedHole(const MeshBase & mesh,
                                              std::set<std::size_t> ids)
  : _center(std::numeric_limits<Real>::max())
{
  // We'll find all the line segments first, then stitch them together
  // afterward
  std::multimap<const Node *, const Node *> hole_edge_map;

  std::vector<boundary_id_type> bcids;

  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      if (elem->dim() == 1)
        {
          if (ids.empty() || ids.count(elem->subdomain_id()))
            {
              hole_edge_map.emplace(elem->node_ptr(0),
                                    elem->node_ptr(1));
              hole_edge_map.emplace(elem->node_ptr(1),
                                    elem->node_ptr(0));
            }
          continue;
        }

      if (elem->dim() == 2)
        {
          const auto ns = elem->n_sides();
          for (auto s : make_range(ns))
            {
              boundary_info.boundary_ids(elem, s, bcids);

              for (auto b : bcids)
                {
                  if (ids.empty() || ids.count(b))
                    {
                      hole_edge_map.emplace(elem->node_ptr(s),
                                            elem->node_ptr((s+1)%ns));
                      hole_edge_map.emplace(elem->node_ptr((s+1)%ns),
                                            elem->node_ptr(s));
                      break;
                    }
                }
            }
        }
    }

  libmesh_error_msg_if
    (hole_edge_map.empty(),
     "No valid hole edges found in mesh!");

  std::vector<const Node *> hole_points
    {hole_edge_map.begin()->first, hole_edge_map.begin()->second};

  // Sort the edges into a connected order
  for (const Node * last = hole_points.front(),
                  *    n = hole_points.back();
                       n != hole_points.front();
                    last = n,
                       n = hole_points.back())
    {
      auto [next_it_begin, next_it_end] = hole_edge_map.equal_range(n);

      libmesh_error_msg_if
        (std::distance(next_it_begin, next_it_end) != 2,
         "Bad edge topology in MeshedHole");

      const Node * next = nullptr;
      for (auto [key, val] : as_range(next_it_begin, next_it_end))
        {
          libmesh_assert_equal_to(key, n);
          libmesh_assert_not_equal_to(val, n);
          if (val == last)
            continue;
          next = val;
        }

      hole_points.push_back(next);
    }

  hole_points.pop_back();

  libmesh_error_msg_if
    (hole_points.size() < 3,
     "Only " << hole_points.size() << " hole edges found in mesh!");

  _points.resize(hole_points.size());
  std::transform(hole_points.begin(),
                 hole_points.end(),
                 _points.begin(),
                 [](const Node * n){ return Point(*n); });

  Real current_area = this->area();

  libmesh_error_msg_if
    (!current_area,
     "Zero-area MeshedHoles are not currently supported");

  // We ordered ourselves backwards?  Whoops, fix that.
  if (current_area < 0)
    std::reverse(_points.begin(), _points.end());
}


unsigned int TriangulatorInterface::MeshedHole::n_points() const
{
  return _points.size();
}


Point TriangulatorInterface::MeshedHole::point(const unsigned int n) const
{
  libmesh_assert_less (n, _points.size());
  return _points[n];
}


Point TriangulatorInterface::MeshedHole::inside() const
{
  // This is expensive to compute, so only do it when we first need it
  if (_center(0) == std::numeric_limits<Real>::max())
    {
      // Start with the vertex average
      _center = std::reduce(_points.begin(), _points.end());

      // Count the number of intersections with a ray to the right,
      // keep track of how far they are
      Point ray_target = _center + Point(1);
      const std::size_t ps = _points.size();
      std::vector<Real> intersection_distances;
      auto find_ray_intersections = [this, ps, &intersection_distances, &ray_target]() {
        for (auto i : make_range(ps))
          {
            const Point & p0 = _points[i],
                        & p1 = _points[(i+1)%ps];
            const Real intersection_distance =
              find_intersection(_center, ray_target, p0, p1);
            if (intersection_distance >= 0)
              intersection_distances.push_back
                (intersection_distance);
          }
      };

      find_ray_intersections();

      // The vertex average isn't on the interior, and we found no
      // intersections to the right?  Try looking to the left.
      if (!intersection_distances.size())
        {
          ray_target = _center - Point(1);
          find_ray_intersections();
        }

      // I'd make this an assert, but I'm not 100% confident we can't
      // get here via some kind of FP error on a weird hole shape.
      libmesh_error_msg_if
        (!intersection_distances.size(),
         "Can't find a center for a MeshedHole!");

      if (!(intersection_distances.size() % 2))
        {
          // Just go from the vertex average to the closest edge
          // intersection, then halfway to the next-closest.

          // Find the nearest first.
          Real min_distance    = std::numeric_limits<Real>::max(),
               second_distance = std::numeric_limits<Real>::max();
          for (Real d : intersection_distances)
            if (d < min_distance)
              {
                second_distance = min_distance;
                min_distance = d;
              }

          const Point ray = ray_target - _center;
          _center += ray * (min_distance + second_distance)/2;
        }
    }
  return _center;
}



} // namespace libMesh
