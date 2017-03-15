// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_BOUNDING_BOX_H
#define LIBMESH_BOUNDING_BOX_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/point.h" // some compilers want the full definition - I think so they can do
// return-value-optimization for BoundingBox'es - BSK

// C++ Includes
#include <vector>
#include <set>
#include <limits>

namespace libMesh
{

/**
 * Defines a Cartesian bounding box by the two
 * corner extremum.
 */
class BoundingBox : public std::pair<Point, Point>
{
public:

  BoundingBox (const Point & new_min,
               const Point & new_max) :
    std::pair<Point, Point>(new_min, new_max)
  {}

  BoundingBox (const std::pair<Point, Point> & bbox) :
    std::pair<Point, Point> (bbox)
  {}

  /**
   * Default constructor sets invalid bounds.
   */
  BoundingBox ()
  {
    this->invalidate();
  }

  /*
   * Sets the bounding box to encompass the universe.
   */
  void invalidate ()
  {
    for (unsigned int i=0; i<LIBMESH_DIM; i++)
      {
        this->first(i)  =  std::numeric_limits<Real>::max();
        this->second(i) = -std::numeric_limits<Real>::max();
      }
  }

  /*
   * Returns a point at the minimum x,y,z coordinates of the box.
   */
  const Point & min() const
  { return this->first; }

  Point & min()
  { return this->first; }

  /*
   * Returns a point at the maximum x,y,z coordinates of the box.
   */
  const Point & max() const
  { return this->second; }

  Point & max()
  { return this->second; }

  /*
   * Returns true iff the other bounding box has a non-empty
   * intersection with this bounding box.
   */
  bool intersects (const BoundingBox &) const;

  bool intersect (const BoundingBox & b) const
  { libmesh_deprecated(); return this->intersects(b); }

  /*
   * Returns true iff the bounding box contains the given point.
   */
  bool contains_point (const Point &) const;

  /*
   * Sets this bounding box to be the intersection with the other
   * bounding box.
   */
  void intersect_with (const BoundingBox &);

  /*
   * Sets this bounding box to be the union with the other
   * bounding box.
   */
  void union_with (const BoundingBox &);

private:
};

} // namespace libMesh


#endif // LIBMESH_BOUNDING_BOX_H
