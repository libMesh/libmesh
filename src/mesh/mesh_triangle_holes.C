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


#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/mesh_triangle_holes.h"

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
// ArbitraryHole member functions
//
TriangulatorInterface::ArbitraryHole::ArbitraryHole(const Point & center,
                                                    const std::vector<Point> & points)
  : _center(center),
    _points(points)
{
  _segment_indices.push_back(0);
  _segment_indices.push_back(points.size());
}

TriangulatorInterface::ArbitraryHole::ArbitraryHole(const Point & center,
                                                    const std::vector<Point> & points,
                                                    const std::vector<unsigned int> & segment_indices)
  : _center(center),
    _points(points),
    _segment_indices(segment_indices)
{}

unsigned int TriangulatorInterface::ArbitraryHole::n_points() const
{
  return _points.size();
}


Point TriangulatorInterface::ArbitraryHole::point(const unsigned int n) const
{
  libmesh_assert_less (n, _points.size());
  return _points[n];
}


Point  TriangulatorInterface::ArbitraryHole::inside() const
{
  return _center;
}

std::vector<unsigned int> TriangulatorInterface::ArbitraryHole::segment_indices() const
{
  return _segment_indices;
}

} // namespace libMesh
