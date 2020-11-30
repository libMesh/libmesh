// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

  /**
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

  /**
   * \returns A point at the minimum x,y,z coordinates of the box.
   */
  const Point & min() const
  { return this->first; }

  Point & min()
  { return this->first; }

  /**
   * \returns A point at the maximum x,y,z coordinates of the box.
   */
  const Point & max() const
  { return this->second; }

  Point & max()
  { return this->second; }

  /**
   * \returns \p true if the other bounding box has a non-empty
   * intersection with this bounding box. Exact floating point <=
   * comparisons are performed.
   */
  bool intersects (const BoundingBox &) const;

  /**
   * \returns \p true if the other bounding box has a non-empty
   * intersection with this bounding box. abstol is an absolute
   * tolerance used to make "fuzzy" comparisons. abstol must be
   * strictly > 0.0, and both BBoxes being compared are "inflated" by
   * abstol in each direction, i.e.
   * (xmin, ymin, zmin) -> (xmin - abstol, ymin - abstol, zmin - abstol)
   * (xmax, ymax, zmax) -> (xmax + abstol, ymax + abstol, zmax + abstol)
   * before the intersection comparisons are made. This approach can
   * be helpful for detecting intersections between two degenerate
   * (planar) bounding boxes that lie in nearly (to within abstol) the
   * same plane and in certain situations should be considered
   * intersecting.
   */
  bool intersects (const BoundingBox &, Real abstol) const;

  /**
   * \returns \p true if the bounding box contains the given point.
   */
  bool contains_point (const Point &) const;

  /**
   * Sets this bounding box to be the intersection with the other
   * bounding box.
   */
  void intersect_with (const BoundingBox &);

  /**
   * Enlarges this bounding box to include the given point
   */
  void union_with (const Point & p);

  /**
   * Sets this bounding box to be the union with the other
   * bounding box.
   */
  void union_with (const BoundingBox &);

  /**
   * Computes the signed distance, d, from a given Point p to this
   * BoundingBox.  The sign convention is:
   * d > 0 if the point is outside the BoundingBox
   * d <= 0 if the point is inside the Bounding Box
   */
  Real signed_distance(const Point & p) const;

  /**
   * Scales each dimension of the bounding box by \p factor.
   *
   * Has no effect for dimensions in which either
   * min(dim) == std::numeric_limits<Real>::max() or
   * max(dim) == -std::numeric_limits<Real>::max(),
   * which is the "invalid" state set by the default
   * constructor and by invalidate().
   */
  void scale(const Real factor);
};



// ------------------------------------------------------------
// BoundingBox class member functions

// BoundingBox::intersects() is about 30% faster when inlined, so its definition
// is here instead of in the source file.
inline
bool
BoundingBox::intersects(const BoundingBox & other_box) const
{
  const libMesh::Point & my_lower = this->first;
  const libMesh::Point & my_upper = this->second;

  const libMesh::Point & other_lower = other_box.first;
  const libMesh::Point & other_upper = other_box.second;

  // Since boxes are tensor products of line intervals it suffices to check
  // that the line segments for each coordinate axis overlap.
  for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
    {
      // Line segments can intersect in two ways:
      // 1. They can overlap.
      // 2. One can be inside the other.
      //
      // In the first case we want to see if either end point of the second
      // line segment lies within the first. In the second case we can simply
      // check that one end point of the first line segment lies in the second
      // line segment. Note that we don't need, in the second case, to do two
      // checks since that case is already covered by the first.
      if (!((my_lower(dir) <= other_lower(dir) &&
             other_lower(dir) <= my_upper(dir)) ||
            (my_lower(dir) <= other_upper(dir) &&
             other_upper(dir) <= my_upper(dir))) &&
          !((other_lower(dir) <= my_lower(dir) &&
             my_lower(dir) <= other_upper(dir))))
        {
          return false;
        }
    }

  return true;
}

inline
bool
BoundingBox::intersects(const BoundingBox & other_box,
                        Real abstol) const
{
  // If you want to use abstol==0, you need to call the "exact"
  // comparison version of the intersects() function.
  libmesh_assert(abstol > 0.);

  BoundingBox expanded_my_box = *this;
  for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
    {
      expanded_my_box.first(dir) -= abstol;
      expanded_my_box.second(dir) += abstol;
    }
  return expanded_my_box.intersects(other_box);
}

inline
void
BoundingBox::union_with(const Point & p)
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      min()(i) = std::min(min()(i), p(i));
      max()(i) = std::max(max()(i), p(i));
    }
}

} // namespace libMesh


#endif // LIBMESH_BOUNDING_BOX_H
