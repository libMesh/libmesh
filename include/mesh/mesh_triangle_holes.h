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

#ifndef LIBMESH_MESH_TRIANGLE_HOLES_H
#define LIBMESH_MESH_TRIANGLE_HOLES_H

#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/point.h"
#include "libmesh/triangulator_interface.h"

// C++ includes

namespace libMesh
{

/**
 * An abstract class for defining a 2-dimensional hole.  We assume that
 * the connectivity of the hole is implicit in the numbering of the points,
 * controlled by segment_indices vector.
 * The size of segment_indices is equal to the number of segments plus one.
 * Each segment has segment_indices[i+1]-segment_indices[i] connected points,
 * with node segment_indices[i] is connected to node segment_indices[i+1],
 * node segment_indices[i+1] is connected to node segment_indices[i+2], etc,
 * and the last node "wraps around" to connect back to node segment_indices[i].
 *
 * \author John W. Peterson
 * \date 2011
 * \brief Class for parameterizing 2D holes to be meshed with Triangle.
 */
class TriangulatorInterface::Hole
{
public:
  /**
   * Constructor
   */
  Hole() = default;

  /**
   * Destructor
   */
  virtual ~Hole() = default;

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

  /**
   * Starting indices of points for a hole with multiple disconnected boundaries.
   */
  virtual std::vector<unsigned int> segment_indices() const
  {
    // default to only one enclosing boundary
    std::vector<unsigned int> seg;
    seg.push_back(0);
    seg.push_back(n_points());
    return seg;
  }
};





/**
 * A concrete instantiation of the Hole class that describes polygonal
 * (triangular, square, pentagonal, ...) holes.
 */
class TriangulatorInterface::PolygonHole : public TriangulatorInterface::Hole
{
public:
  /**
   * Constructor specifying the center, radius, and number of
   * points which comprise the hole.  The points will all lie
   * on a circle of radius r.
   */
  PolygonHole(const Point & center, Real radius, unsigned int n_points);

  virtual unsigned int n_points() const override;

  virtual Point point(const unsigned int n) const override;

  virtual Point inside() const override;

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
class TriangulatorInterface::ArbitraryHole : public TriangulatorInterface::Hole
{
public:
  /**
   * The constructor requires a point which lies in the interior of the hole
   * and a reference to a vector of Points defining the hole.
   */
  ArbitraryHole(const Point & center,
                const std::vector<Point> & points);

  ArbitraryHole(const Point & center,
                const std::vector<Point> & points,
                const std::vector<unsigned int> & segment_indices);

  virtual unsigned int n_points() const override;

  virtual Point point(const unsigned int n) const override;

  virtual Point inside() const override;

  virtual std::vector<unsigned int> segment_indices() const override;

private:
  /**
   * arbitrary (x,y) location inside the hole
   */
  const Point _center;

  /**
   * Reference to the vector of points which makes up
   * the hole.
   */
  const std::vector<Point> _points;

  std::vector<unsigned int> _segment_indices;
};





/**
 * A class for defining a 2-dimensional region for Triangle.
 */
class TriangulatorInterface::Region
{
public:
  /**
   * Constructor
   */
  Region(const Point & center,
         Real attribute = 0,
         Real max_area = -std::numeric_limits<Real>::max())
  : _center(center),
    _attribute(attribute),
    _max_area(max_area) {}

  Point inside() const { return _center; }

  Real attribute() const { return _attribute; }

  Real max_area() const { return _max_area; }

private:
  /**
   * Arbitrary (x,y) location inside the region
   */
  const Point _center;

  /**
   * Attribute of the region
   * Default value for attribute is zero.
   */
  Real _attribute;

  /**
   * Maximum area that is allowed for all triangles in this region
   * Default negative maximum area means that no area constraint will be imposed.
   */
  Real _max_area;
};

} // namespace libMesh

#endif // LIBMESH_MESH_TRIANGLE_HOLES_H
