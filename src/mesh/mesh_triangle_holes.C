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


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local includes
#include "libmesh/mesh_triangle_holes.h"

namespace libMesh
{

//
// PolygonHole member functions
//
TriangleInterface::PolygonHole::PolygonHole(const Point & center,
                                            Real radius,
                                            unsigned int n_points_in) :
  _center(center),
  _radius(radius),
  _n_points(n_points_in)
{}


unsigned int TriangleInterface::PolygonHole::n_points() const
{
  return _n_points;
}

Point TriangleInterface::PolygonHole::point(const unsigned int n) const
{
  // The nth point lies at the angle theta = 2 * pi * n / _n_points
  const Real theta = static_cast<Real>(n) * 2.0 * libMesh::pi / static_cast<Real>(_n_points);

  return Point(_center(0) + _radius*std::cos(theta), // x=r*cos(theta)
               _center(1) + _radius*std::sin(theta), // y=r*sin(theta)
               0.);
}



Point TriangleInterface::PolygonHole::inside() const
{
  // The center of the hole is definitely inside.
  return _center;
}



//
// ArbitraryHole member functions
//
TriangleInterface::ArbitraryHole::ArbitraryHole(const Point & center,
                                                const std::vector<Point> & points)
  : _center(center),
    _points(points)
{}


unsigned int TriangleInterface::ArbitraryHole::n_points() const
{
  return _points.size();
}


Point TriangleInterface::ArbitraryHole::point(const unsigned int n) const
{
  libmesh_assert_less (n, _points.size());
  return _points[n];
}


Point  TriangleInterface::ArbitraryHole::inside() const
{
  return _center;
}


} // namespace libMesh


#endif // LIBMESH_HAVE_TRIANGLE
