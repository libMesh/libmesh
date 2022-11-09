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



#ifndef LIBMESH_INTERSECTIONTOOLS_H
#define LIBMESH_INTERSECTIONTOOLS_H

#include "libmesh/libmesh_common.h"

namespace libMesh
{

class Point;

namespace IntersectionTools
{

/**
 * Enum that represents the result of the \p within_segment() methods.
*/
enum WithinSegmentResult : int {
                         NOT_WITHIN = 0,
                         AT_BEGINNING = 1,
                         AT_END = 2,
                         BETWEEN = 3 };

/**
 * Checks whether or not a point is within a line segment
 * @param s1 The first point on the segment
 * @param s2 The second point on the segment
 * @param length The segment length (for optimization if it's already computed)
 * @param p The point
 * @param tol The tolerance to use
 * @return Enum denoting whether or not the point is not within, at s1,
 * at s2, or between [s1, s2]
 */
WithinSegmentResult within_segment(const Point & s1,
                                   const Point & s2,
                                   const Real length,
                                   const Point & p,
                                   const Real tol = TOLERANCE);

/**
 * Checks whether or not a point is within a line segment
 * @param s1 The first point on the segment
 * @param s2 The second point on the segment
 * @param p The point
 * @param tol The tolerance to use
 * @return Enum denoting whether or not the point is not within, at s1,
 * at s2, or between [s1, s2]
 */
WithinSegmentResult within_segment(const Point & s1,
                                   const Point & s2,
                                   const Point & p,
                                   const Real tol = TOLERANCE);

/**
 * Checks whether or not the given points are collinear
 * @param p1 The first point
 * @param p2 The second point
 * @param p3 The third point
 * @param tol The tolerance to use
 * @return Whether or not the given points are collinear
 */
bool collinear(const Point & p1,
               const Point & p2,
               const Point & p3,
               const Real tol = TOLERANCE);

} // namespace IntersectionTools
} // namespace libMesh

#endif // LIBMESH_INTERSECTIONTOOLS_H
