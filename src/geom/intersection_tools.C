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



#include "libmesh/intersection_tools.h"

#include "libmesh/point.h"
#include "libmesh/int_range.h"

namespace libMesh::IntersectionTools
{

WithinSegmentResult within_segment(const Point & s1,
                                   const Point & s2,
                                   const Real length,
                                   const Point & p,
                                   const Real tol)
{
  libmesh_assert(!s1.absolute_fuzzy_equals(s2, tol));
  libmesh_assert_less(std::abs((s1 - s2).norm() - length), tol);

  const auto diff1 = p - s1;
  const auto diff2 = p - s2;
  const auto tol_scaled = tol * length;

  if (diff1 * diff2 > tol_scaled)
    return NOT_WITHIN;

  const auto diff1_norm = diff1.norm();
  if (diff1_norm < tol_scaled)
    return AT_BEGINNING;

  const auto diff2_norm = diff2.norm();
  if (diff2_norm < tol_scaled)
    return AT_END;

  // whether or not p is _between_ [s1, s2]
  if (std::abs(diff1_norm + diff2_norm - length) < tol_scaled)
    return BETWEEN;
  return NOT_WITHIN;
}

WithinSegmentResult within_segment(const Point & s1,
                                   const Point & s2,
                                   const Point & p,
                                   const Real tol)
{
  return within_segment(s1, s2, (s1 - s2).norm(), p, tol);
}

bool collinear(const Point & p1,
               const Point & p2,
               const Point & p3,
               const Real tol)
{
  // (p1 -> p2) X (p1 - > p3) should be == the zero vector
  const auto p1p2_cross_p1p3 = (p2 - p1).cross(p3 - p1);
  for (const auto i : make_range(LIBMESH_DIM))
    if (std::abs(p1p2_cross_p1p3(i)) > tol)
      return false;
  return true;
}

} // namespace libMesh::IntersectionTools


