// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_MESH_TRIANGLE_HOLES_H
#define LIBMESH_MESH_TRIANGLE_HOLES_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local includes
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/point.h"

// C++ includes

namespace libMesh
{

/**
 * An abstract class for defining a 2-dimensional hole.  We assume that
 * the connectivity of the hole is implicit in the numbering of the points,
 * i.e. node 0 is connected to node 1, node 1 is connected to node 2, etc,
 * and the last node "wraps around" to connect back to node 0.
 */
class TriangleInterface::Hole
{
public:
  /**
   * Constructor
   */
  Hole() {}

  /**
   * Destructor
   */
  virtual ~Hole() {}

  /**
   * The number of geometric points which define the hole.
   */
  virtual unsigned int n_points() const = 0;

  /**
   * Return the nth point defining the hole.
   */
  virtual Point point(const unsigned int n) const = 0;

  /**
   * Return an (arbitrary) point which lies inside the hole.
   */
  virtual Point inside() const = 0;
};





/**
 * A concrete instantiation of the Hole class that describes polygonal
 * (triangular, square, pentagonal, ...) holes.
 */
class TriangleInterface::PolygonHole : public TriangleInterface::Hole
{
public:
  /**
   * Constructor specifying the center, radius, and number of
   * points which comprise the hole.  The points will all lie
   * on a circle of radius r.
   */
  PolygonHole(const Point & center, Real radius, unsigned int n_points);

  virtual unsigned int n_points() const libmesh_override;

  virtual Point point(const unsigned int n) const libmesh_override;

  virtual Point inside() const libmesh_override;

private:
  /**
   * (x,y) location of the center of the hole
   */
  Point _center;

  /**
   * circular hole radius
   */
  Real _radius;

  /**
   * number of points used to describe the hole.  The actual
   * points can be generated knowing the center and radius.
   * For example, n_points=3 would generate a triangular hole.
   */
  unsigned int _n_points;
};




/**
 * Another concrete instantiation of the hole, this one should
 * be sufficiently general for most non-polygonal purposes.  The user
 * supplies, at the time of construction, a reference to a vector
 * of Points which defines the hole (in order of connectivity) and
 * an arbitrary Point which lies inside the hole.
 */
class TriangleInterface::ArbitraryHole : public TriangleInterface::Hole
{
public:
  /**
   * The constructor requires a point which lies in the interior of the hole
   * and a reference to a vector of Points defining the hole.
   */
  ArbitraryHole(const Point & center,
                const std::vector<Point> & points);

  virtual unsigned int n_points() const libmesh_override;

  virtual Point point(const unsigned int n) const libmesh_override;

  virtual Point inside() const libmesh_override;

private:
  /**
   * arbitrary (x,y) location inside the hole
   */
  Point _center;

  /**
   * Reference to the vector of points which makes up
   * the hole.
   */
  const std::vector<Point> & _points;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE

#endif // LIBMESH_MESH_TRIANGLE_HOLES_H
