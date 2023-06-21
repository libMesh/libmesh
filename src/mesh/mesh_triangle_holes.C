// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/mesh_serializer.h"
#include "libmesh/node.h"
#include "libmesh/parallel_algebra.h" // Packing<Point>
#include "libmesh/simple_range.h"

// TIMPI includes
#include "timpi/parallel_implementation.h" // broadcast

// C++ includes
#include <algorithm>

namespace
{
  using namespace libMesh;

  int signof(Real val) {
    return (0 < val) - (val < 0);
  }

  // Return 1 iff counter-clockwise turn
  // Return -1 iff clockwise turn
  // Return 0 iff collinear
  int orientation(const Point & p0,
                  const Point & p1,
                  const Point & p2)
  {
    const double detleft  = (p0(0)-p2(0))*(p1(1)-p2(1));
    const double detright = (p0(1)-p2(1))*(p1(0)-p2(0));

    return signof(detleft - detright);
  }

  // Same, but for the ray target as it goes to infinity
  int ray_orientation(const Point & p0,
                      const Point & p1,
                      const Point & source,
                      const Point & ray_target)
  {
    const Point rayvec = ray_target - source;
    const Point edgevec = p1 - p0;
    const double det = edgevec(0)*rayvec(1)-edgevec(1)*rayvec(0);

    return signof(det);
  }

  bool is_intersection(const Point & source,
                       const Point & ray_target,
                       const Point & edge_pt0,
                       const Point & edge_pt1)
  {
    int orient_st0 = orientation(source, ray_target, edge_pt0);
    int orient_st1 = orientation(source, ray_target, edge_pt1);
    int orient_edge_s = orientation(edge_pt0, edge_pt1, source);
    int orient_edge_t = ray_orientation(edge_pt0, edge_pt1, source, ray_target);

    // Intersection on interior
    if ((orient_st0 == -orient_st1) &&
        (orient_edge_s != orient_edge_t))
      return true;

    // Ray intersects edge_pt1
    if (orient_st1 == 0)
      return true;

    // Source is on line; we don't count that
    // if (orient_edge_s == 0)
    // Ray is parallel to edge; no intersection;
    // if (orient_edge_t == 0)
    // Ray intersects edge_pt0; we don't count that
    // if (orient_st0 == 0)

    return false;
  }

  // Returns a positive distance iff the ray from source in the
  // direction of ray_target intersects the edge from pt0
  // (non-inclusive) to pt1 (inclusive), -1 otherwise.
  //
  // If the intersection is a "glancing" one at a corner, return -1.
  Real find_intersection(const Point & source,
                         const Point & ray_target_0,
                         const Point & edge_pt0,
                         const Point & edge_pt1,
                         const Point & edge_pt2)
  {
    // Add a random small shift to the source and three points
    const Point ray_target = ray_target_0 + Point((Real)(rand() % 5000 + 5000) / 5000.0 * libMesh::TOLERANCE * libMesh::TOLERANCE, (Real)(rand() % 5000 + 5000) / 5000.0 * libMesh::TOLERANCE * libMesh::TOLERANCE, 0.0);

    // Quick and more numerically stable check
    if (!is_intersection(source, ray_target, edge_pt0, edge_pt1))
      return -1;

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
        // There's an intersection right on a vertex?  We'll count it
        // if and only if it isn't a "double-intersection", if the
        // *next* edge in line is on the other side of our ray.
        if (!t)
          {
            const Real prevdx = edge_pt0(0)-ray_target(0),
                       prevdy = edge_pt0(1)-ray_target(1);
            const Real p_num = prevdx * raydy -
                               prevdy * raydx;

            const Real nextdx = edge_pt2(0)-ray_target(0),
                       nextdy = edge_pt2(1)-ray_target(1);
            const Real n_num = nextdx * raydy -
                               nextdy * raydx;

            if (signof(p_num) != -signof(n_num))
              return -1;
          }

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
  return this->areavec().norm() / 2;
}


RealGradient TriangulatorInterface::Hole::areavec() const
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
  // regardless of which plane they're in.  If you're in the XY plane,
  // then the standard counter-clockwise hole point ordering gives you
  // a positive areavec(2);

  RealGradient areavec = 0;

  for (unsigned int i=2; i != np; ++i)
    {
      const Point e_0im = this->point(i-1) - p0,
                  e_0i  = this->point(i) - p0;

      areavec += e_0i.cross(e_0im);
    }

  return areavec;
}



std::vector<Real>
TriangulatorInterface::Hole::find_ray_intersections(Point ray_start,
                                                    Point ray_target) const
{
  const auto np = this->n_points();

  std::vector<Real> intersection_distances;

  for (auto i : make_range(np))
    {
      const Point & p0 = this->point(i),
                  & p1 = this->point((i+1)%np),
                  & p2 = this->point((i+2)%np);
      const Real intersection_distance =
        find_intersection(ray_start, ray_target, p0, p1, p2);
      if (intersection_distance >= 0)
        intersection_distances.push_back
          (intersection_distance);
    }

  return intersection_distances;
}



Point TriangulatorInterface::Hole::calculate_inside_point() const
{
  // Start with the vertex average

  // Turns out "I'm a fully compliant C++17 compiler!" doesn't
  // mean "I have a full C++17 standard library!"
  // inside = std::reduce(points.begin(), points.end());
  Point inside = 0;
  for (auto i : make_range(this->n_points()))
    inside += this->point(i);

  inside /= this->n_points();

  // Count the number of intersections with a ray to the right,
  // keep track of how far they are
  Point ray_target = inside + Point(1);
  std::vector<Real> intersection_distances =
    this->find_ray_intersections(inside, ray_target);

  // The vertex average isn't on the interior, and we found no
  // intersections to the right?  Try looking to the left.
  if (!intersection_distances.size())
    {
      ray_target = inside - Point(1);
      intersection_distances =
        this->find_ray_intersections(inside, ray_target);
    }

  // I'd make this an assert, but I'm not 100% confident we can't
  // get here via some kind of FP error on a weird hole shape.
  libmesh_error_msg_if
    (!intersection_distances.size(),
     "Can't find a center for a MeshedHole!");

  if (intersection_distances.size() % 2)
    return inside;

  // The vertex average is outside.  So go from the vertex average to
  // the closest edge intersection, then halfway to the next-closest.

  // Find the nearest first.
  Real min_distance    = std::numeric_limits<Real>::max(),
       second_distance = std::numeric_limits<Real>::max();
  for (Real d : intersection_distances)
    if (d < min_distance)
      {
        second_distance = min_distance;
        min_distance = d;
      }

  const Point ray = ray_target - inside;
  inside += ray * (min_distance + second_distance)/2;

  return inside;
}


bool TriangulatorInterface::Hole::contains(Point p) const
{
  // Count the number of intersections with a ray to the right,
  // keep track of how far they are
  Point ray_target = p + Point(1);
  std::vector<Real> intersection_distances =
    this->find_ray_intersections(p, ray_target);

  // Odd number of intersections == we're inside
  // Even number == we're outside
  return intersection_distances.size() % 2;
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


TriangulatorInterface::ArbitraryHole::ArbitraryHole(std::vector<Point> points)
  : _points(std::move(points))
{
  _segment_indices.push_back(0);
  _segment_indices.push_back(_points.size());
  _center = this->calculate_inside_point();
}


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
  // We'll want to do this on one processor and broadcast to the rest;
  // otherwise we can get out of sync by doing things like using
  // pointers as keys.
  libmesh_parallel_only(mesh.comm());

  MeshSerializer serial(const_cast<MeshBase &>(mesh),
                        /* serial */ true, /* only proc 0 */ true);

  // Try to keep in sync even if we throw an error on proc 0, so we
  // can examine errors in our unit tests in parallel too.
  std::string error_reported;

  auto report_error = [&mesh, &error_reported](std::string er) {
    error_reported = std::move(er);
    mesh.comm().broadcast(error_reported);
    libmesh_error_msg(error_reported);
  };

  if (mesh.processor_id() != 0)
    {
      // Make sure proc 0 didn't just fail
      mesh.comm().broadcast(error_reported);
      libmesh_error_msg_if(!error_reported.empty(), error_reported);

      // Receive the points proc 0 will send later
      mesh.comm().broadcast(_points);
      return;
    }

  // We'll find all the line segments first, then stitch them together
  // afterward.  If the line segments come from 2D element sides then
  // we'll label their edge_type as "1" for clockwise orientation
  // around the element or "2" for CCW, to make it easier to detect
  // and scream about cases where we have a disconnected outer
  // boundary.
  std::multimap<const Node *, std::pair<const Node *, int>> hole_edge_map;

  std::vector<boundary_id_type> bcids;

  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      if (elem->dim() == 1)
        {
          if (ids.empty() || ids.count(elem->subdomain_id()))
            {
              hole_edge_map.emplace(elem->node_ptr(0),
                                    std::make_pair(elem->node_ptr(1),
                                                   /*edge*/ 0));
              hole_edge_map.emplace(elem->node_ptr(1),
                                    std::make_pair(elem->node_ptr(0),
                                                   /*edge*/ 0));
            }
          continue;
        }

      if (elem->dim() == 2)
        {
          const auto ns = elem->n_sides();
          for (auto s : make_range(ns))
            {
              boundary_info.boundary_ids(elem, s, bcids);

              bool add_edge = false;
              if (!elem->neighbor_ptr(s) && ids.empty())
                add_edge = true;

              if (!add_edge)
                for (auto b : bcids)
                  if (ids.count(b))
                    add_edge = true;

              if (add_edge)
                {
                  hole_edge_map.emplace(elem->node_ptr(s),
                                        std::make_pair(elem->node_ptr((s+1)%ns),
                                                       /*counter-CW*/ 2));
                  // Do we really need to support flipped 2D elements?
                  hole_edge_map.emplace(elem->node_ptr((s+1)%ns),
                                        std::make_pair(elem->node_ptr(s),
                                                       /*clockwise*/ 1));
                  continue;
                }
            }
        }
    }

  if (hole_edge_map.empty())
    report_error("No valid hole edges found in mesh!");

  // Function to pull a vector of points out of the map; a loop of
  // edges connecting these points defines a hole boundary.  If the
  // mesh has multiple boundaries (e.g. because it had holes itself),
  // then a random vector will be extracted; this function will be
  // called multiple times so that the various options can be
  // compared.  We choose the largest option.
  auto extract_edge_vector = [&report_error, &hole_edge_map]() {
    // Start with any edge
    std::pair<std::vector<const Node *>, int> hole_points_and_edge_type
    {{hole_edge_map.begin()->first, hole_edge_map.begin()->second.first},
     hole_edge_map.begin()->second.second};

    int & edge_type = hole_points_and_edge_type.second;
    auto & hole_points = hole_points_and_edge_type.first;

    // We won't be needing to search for this edge
    hole_edge_map.erase(hole_points.front());

    // Sort the remaining edges into a connected order
    for (const Node * last = hole_points.front(),
                    *    n = hole_points.back();
                         n != hole_points.front();
                      last = n,
                         n = hole_points.back())
      {
        auto [next_it_begin, next_it_end] = hole_edge_map.equal_range(n);

        if (std::distance(next_it_begin, next_it_end) != 2)
          report_error("Bad edge topology found by MeshedHole");

        const Node * next = nullptr;
        for (const auto & [key, val] : as_range(next_it_begin, next_it_end))
          {
            libmesh_assert_equal_to(key, n);
            libmesh_ignore(key);
            libmesh_assert_not_equal_to(val.first, n);

            // Don't go backwards on the edge we just traversed
            if (val.first == last)
              continue;

            // We can support mixes of Edge and Tri-side edges, but we
            // can't do proper error detection on flipped triangles.
            if (val.second != edge_type &&
                val.second != 0)
              {
                if (!edge_type)
                  edge_type = val.second;
                else
                  report_error("MeshedHole sees inconsistent triangle orientations on boundary");
              }
            next = val.first;
          }

        // We should never hit the same n twice!
        hole_edge_map.erase(next_it_begin, next_it_end);

        hole_points.push_back(next);
      }

    hole_points.pop_back();

    return hole_points_and_edge_type;
  };

  /*
   * If it's not obvious which loop we find is really the loop we
   * want, then we should die with a nice error message.
   */
  int n_negative_areas = 0,
      n_positive_areas = 0,
      n_edgeelem_loops = 0;

  std::vector<const Node *> outer_hole_points;
  int outer_edge_type = -1;
  Real twice_outer_area = 0,
       abs_twice_outer_area = 0;

#ifdef DEBUG
  // Area and edge type, for error reporting
  std::vector<std::pair<Real, int>> areas;
#endif

  while (!hole_edge_map.empty()) {
    auto [hole_points, edge_type] = extract_edge_vector();

    if (edge_type == 0)
    {
      ++n_edgeelem_loops;
      if (n_edgeelem_loops > 1)
        report_error("MeshedHole is confused by multiple loops of Edge elements");
      if (n_positive_areas || n_negative_areas)
        report_error("MeshedHole is confused by meshes with both Edge and 2D-side boundaries");
    }

    const std::size_t n_hole_points = hole_points.size();
    if (n_hole_points < 3)
      report_error("Loop with only " + std::to_string(n_hole_points) +
                   " hole edges found in mesh!");

    Real twice_this_area = 0;
    const Point p0 = *hole_points[0];
    for (unsigned int i=2; i != n_hole_points; ++i)
      {
        const Point e_0im = *hole_points[i-1] - p0,
                    e_0i  = *hole_points[i] - p0;

        twice_this_area += e_0i.cross(e_0im)(2);
      }

    auto abs_twice_this_area = std::abs(twice_this_area);

    if (((twice_this_area > 0) && edge_type == 2) ||
        ((twice_this_area < 0) && edge_type == 1))
      ++n_positive_areas;
    else
      ++n_negative_areas;

#ifdef DEBUG
    areas.push_back({twice_this_area/2,edge_type});
#endif

    if (abs_twice_this_area > abs_twice_outer_area)
      {
        twice_outer_area = twice_this_area;
        abs_twice_outer_area = abs_twice_this_area;
        outer_hole_points = std::move(hole_points);
        outer_edge_type = edge_type;
      }
  }

  _points.resize(outer_hole_points.size());
  std::transform(outer_hole_points.begin(),
                 outer_hole_points.end(),
                 _points.begin(),
                 [](const Node * n){ return Point(*n); });

  if (!twice_outer_area)
    report_error("Zero-area MeshedHoles are not currently supported");

  // We ordered ourselves counter-clockwise?  But a hole is expected
  // to be clockwise, so use the reverse order.
  if (twice_outer_area > 0)
    std::reverse(_points.begin(), _points.end());

#ifdef DEBUG
  auto print_areas = [areas](){
    libMesh::out << "Found boundary areas:\n";
    static const std::vector<std::string> edgenames {"E","CW","CCW"};
    for (auto area : areas)
      libMesh::out << '(' << edgenames[area.second] << ' ' <<
        area.first << ')';
    libMesh::out << std::endl;
  };
#else
  auto print_areas = [](){};
#endif

  if (((twice_outer_area > 0) && outer_edge_type == 2) ||
      ((twice_outer_area < 0) && outer_edge_type == 1))
    {
      if (n_positive_areas > 1)
        {
          print_areas();
          report_error("MeshedHole found " +
                       std::to_string(n_positive_areas) +
                       " counter-clockwise boundaries and cannot choose one!");
        }

    }
  else if (outer_edge_type != 0)
    {
      if (n_negative_areas > 1)
        {
          print_areas();
          report_error("MeshedHole found " +
                       std::to_string(n_positive_areas) +
                       " clockwise boundaries and cannot choose one!");
        }

    }

  // Hey, no errors!  Broadcast that empty string.
  mesh.comm().broadcast(error_reported);
  mesh.comm().broadcast(_points);
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
    _center = this->calculate_inside_point();

  return _center;
}



} // namespace libMesh
