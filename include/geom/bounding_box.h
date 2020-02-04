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
#include <algorithm>
#include <array>

namespace libMesh
{

/**
 * Defines a Cartesian bounding box by the two
 * corner extremum.
 */
template <typename RealType = Real>
class BoundingBoxTempl : public std::pair<PointTempl<RealType>, PointTempl<RealType>>
{
public:
  typedef BoundingBoxTempl<RealType> BoundingBox;
  typedef PointTempl<RealType> Point;

  BoundingBoxTempl (const Point & new_min,
                    const Point & new_max) :
    std::pair<Point, Point>(new_min, new_max)
  {}

  BoundingBoxTempl (const std::pair<Point, Point> & bbox) :
    std::pair<Point, Point> (bbox)
  {}

  /**
   * Default constructor sets invalid bounds.
   */
  BoundingBoxTempl ()
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
   * \returns \p true if the other bounding box has a non-empty
   * intersection with this bounding box.
   *
   * \deprecated Use the BoundingBox::intersects() function instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  bool intersect (const BoundingBox & b) const
  { libmesh_deprecated(); return this->intersects(b); }
#endif

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
  RealType signed_distance(const Point & p) const;
};



// ------------------------------------------------------------
// BoundingBox class member functions
template <typename RealType>
inline
void
BoundingBoxTempl<RealType>::union_with(const Point & p)
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      min()(i) = std::min(min()(i), p(i));
      max()(i) = std::max(max()(i), p(i));
    }
}

namespace {
// Small helper function to make intersects() more readable.
template <typename RealType>
bool is_between(const RealType & min, const RealType & check, const RealType & max)
{
  return min <= check && check <= max;
}
}

template <typename RealType>
bool BoundingBoxTempl<RealType>::intersects (const BoundingBox & other_box) const
{
  // Make local variables first to make things more clear in a moment
  const auto & my_min_x = this->first(0);
  const auto & my_max_x = this->second(0);
  const auto & other_min_x = other_box.first(0);
  const auto & other_max_x = other_box.second(0);

  const bool x_int = is_between(my_min_x, other_min_x, my_max_x) || is_between(my_min_x, other_max_x, my_max_x) ||
    is_between(other_min_x, my_min_x, other_max_x) || is_between(other_min_x, my_max_x, other_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  const auto & my_min_y = this->first(1);
  const auto & my_max_y = this->second(1);
  const auto & other_min_y = other_box.first(1);
  const auto & other_max_y = other_box.second(1);

  const bool y_int = is_between(my_min_y, other_min_y, my_max_y) || is_between(my_min_y, other_max_y, my_max_y) ||
    is_between(other_min_y, my_min_y, other_max_y) || is_between(other_min_y, my_max_y, other_max_y);

  intersection_true = intersection_true && y_int;
#endif

#if LIBMESH_DIM > 2
  const auto & my_min_z = this->first(2);
  const auto & my_max_z = this->second(2);
  const auto & other_min_z = other_box.first(2);
  const auto & other_max_z = other_box.second(2);

  const bool z_int = is_between(my_min_z, other_min_z, my_max_z) || is_between(my_min_z, other_max_z, my_max_z) ||
    is_between(other_min_z, my_min_z, other_max_z) || is_between(other_min_z, my_max_z, other_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

template <typename RealType>
bool BoundingBoxTempl<RealType>::intersects (const BoundingBox & other_box,
                                             Real abstol) const
{
  // If you want to use abstol==0, you need to call the "exact"
  // comparison version of the intersects() function.
  libmesh_assert(abstol > 0.);

  // Make local variables first to make things more clear in a moment
  const auto & my_min_x = this->first(0);
  const auto & my_max_x = this->second(0);
  const auto & ot_min_x = other_box.first(0);
  const auto & ot_max_x = other_box.second(0);

  const bool x_int =
    is_between(my_min_x - abstol, ot_min_x, my_max_x + abstol) ||
    is_between(my_min_x - abstol, ot_max_x, my_max_x + abstol) ||
    is_between(ot_min_x - abstol, my_min_x, ot_max_x + abstol) ||
    is_between(ot_min_x - abstol, my_max_x, ot_max_x + abstol);

  bool intersection_true = x_int;

  if (!intersection_true)
    return false;

#if LIBMESH_DIM > 1
  const auto & my_min_y = this->first(1);
  const auto & my_max_y = this->second(1);
  const auto & ot_min_y = other_box.first(1);
  const auto & ot_max_y = other_box.second(1);

  const bool y_int =
    is_between(my_min_y - abstol, ot_min_y, my_max_y + abstol) ||
    is_between(my_min_y - abstol, ot_max_y, my_max_y + abstol) ||
    is_between(ot_min_y - abstol, my_min_y, ot_max_y + abstol) ||
    is_between(ot_min_y - abstol, my_max_y, ot_max_y + abstol);

  intersection_true = intersection_true && y_int;

  if (!intersection_true)
    return false;
#endif

#if LIBMESH_DIM > 2
  const auto & my_min_z = this->first(2);
  const auto & my_max_z = this->second(2);
  const auto & ot_min_z = other_box.first(2);
  const auto & ot_max_z = other_box.second(2);

  const bool z_int =
    is_between(my_min_z - abstol, ot_min_z, my_max_z + abstol) ||
    is_between(my_min_z - abstol, ot_max_z, my_max_z + abstol) ||
    is_between(ot_min_z - abstol, my_min_z, ot_max_z + abstol) ||
    is_between(ot_min_z - abstol, my_max_z, ot_max_z + abstol);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

template <typename RealType>
bool BoundingBoxTempl<RealType>::contains_point (const Point & p) const
{
  // Make local variables first to make things more clear in a moment
  const auto & my_min_x = this->first(0);
  const auto & my_max_x = this->second(0);
  bool x_int = is_between(my_min_x, p(0), my_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  const auto & my_min_y = this->first(1);
  const auto & my_max_y = this->second(1);
  bool y_int = is_between(my_min_y, p(1), my_max_y);

  intersection_true = intersection_true && y_int;
#endif


#if LIBMESH_DIM > 2
  const auto & my_min_z = this->first(2);
  const auto & my_max_z = this->second(2);
  bool z_int = is_between(my_min_z, p(2), my_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}


template <typename RealType>
void BoundingBoxTempl<RealType>::intersect_with (const BoundingBox & other_box)
{
  this->first(0)  = std::max(this->first(0),  other_box.first(0));
  this->second(0) = std::min(this->second(0), other_box.second(0));

#if LIBMESH_DIM > 1
  this->first(1)  = std::max(this->first(1),  other_box.first(1));
  this->second(1) = std::min(this->second(1), other_box.second(1));
#endif

#if LIBMESH_DIM > 2
  this->first(2)  = std::max(this->first(2),  other_box.first(2));
  this->second(2) = std::min(this->second(2), other_box.second(2));
#endif
}


template <typename RealType>
void BoundingBoxTempl<RealType>::union_with (const BoundingBox & other_box)
{
  this->first(0)  = std::min(this->first(0),  other_box.first(0));
  this->second(0) = std::max(this->second(0), other_box.second(0));

#if LIBMESH_DIM > 1
  this->first(1)  = std::min(this->first(1),  other_box.first(1));
  this->second(1) = std::max(this->second(1), other_box.second(1));
#endif

#if LIBMESH_DIM > 2
  this->first(2)  = std::min(this->first(2),  other_box.first(2));
  this->second(2) = std::max(this->second(2), other_box.second(2));
#endif
}



template <typename RealType>
RealType BoundingBoxTempl<RealType>::signed_distance(const Point & p) const
{
  if (contains_point(p))
    {
      // Sign convention: if Point is inside the bbox, the distance is
      // negative. We then find the smallest distance to the different
      // sides of the box and return that.
      RealType min_dist = std::numeric_limits<Real>::max();

      for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
        {
          min_dist = std::min(min_dist, std::abs(p(dir) - this->second(dir)));
          min_dist = std::min(min_dist, std::abs(p(dir) - this->first(dir)));
        }

      return -min_dist;
    }
  else // p is outside the box
    {
      RealType dx[3] = {0., 0., 0.};

      // Compute distance "above"/"below" the box in each
      // direction. If the point is somewhere in between the (min,
      // max) values of the box, dx is 0.
      for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
        {
          if (p(dir) > this->second(dir))
            dx[dir] = p(dir) - this->second(dir);
          else if (p(dir) < this->first(dir))
            dx[dir] = p(dir) - this->first(dir);
        }

      return std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }
}

typedef BoundingBoxTempl<Real> BoundingBox;

} // namespace libMesh


#endif // LIBMESH_BOUNDING_BOX_H
