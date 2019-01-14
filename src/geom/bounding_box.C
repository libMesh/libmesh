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



// Local includes
#include "libmesh/bounding_box.h"

// C++ includes
#include <algorithm> // std::min_element
#include <array>

namespace libMesh
{
// Small helper function to make intersects() more readable.
bool is_between(Real min, Real check, Real max)
{
  return min <= check && check <= max;
}

bool BoundingBox::intersects (const BoundingBox & other_box) const
{
  // Make local variables first to make things more clear in a moment
  const Real & my_min_x = this->first(0);
  const Real & my_max_x = this->second(0);
  const Real & other_min_x = other_box.first(0);
  const Real & other_max_x = other_box.second(0);

  const bool x_int = is_between(my_min_x, other_min_x, my_max_x) || is_between(my_min_x, other_max_x, my_max_x) ||
    is_between(other_min_x, my_min_x, other_max_x) || is_between(other_min_x, my_max_x, other_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  const Real & my_min_y = this->first(1);
  const Real & my_max_y = this->second(1);
  const Real & other_min_y = other_box.first(1);
  const Real & other_max_y = other_box.second(1);

  const bool y_int = is_between(my_min_y, other_min_y, my_max_y) || is_between(my_min_y, other_max_y, my_max_y) ||
    is_between(other_min_y, my_min_y, other_max_y) || is_between(other_min_y, my_max_y, other_max_y);

  intersection_true = intersection_true && y_int;
#endif

#if LIBMESH_DIM > 2
  const Real & my_min_z = this->first(2);
  const Real & my_max_z = this->second(2);
  const Real & other_min_z = other_box.first(2);
  const Real & other_max_z = other_box.second(2);

  const bool z_int = is_between(my_min_z, other_min_z, my_max_z) || is_between(my_min_z, other_max_z, my_max_z) ||
    is_between(other_min_z, my_min_z, other_max_z) || is_between(other_min_z, my_max_z, other_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

bool BoundingBox::intersects (const BoundingBox & other_box,
                              Real abstol) const
{
  // If you want to use abstol==0, you need to call the "exact"
  // comparison version of the intersects() function.
  libmesh_assert(abstol > 0.);

  // Make local variables first to make things more clear in a moment
  const Real & my_min_x = this->first(0);
  const Real & my_max_x = this->second(0);
  const Real & ot_min_x = other_box.first(0);
  const Real & ot_max_x = other_box.second(0);

  const bool x_int =
    is_between(my_min_x - abstol, ot_min_x, my_max_x + abstol) ||
    is_between(my_min_x - abstol, ot_max_x, my_max_x + abstol) ||
    is_between(ot_min_x - abstol, my_min_x, ot_max_x + abstol) ||
    is_between(ot_min_x - abstol, my_max_x, ot_max_x + abstol);

  bool intersection_true = x_int;

  if (!intersection_true)
    return false;

#if LIBMESH_DIM > 1
  const Real & my_min_y = this->first(1);
  const Real & my_max_y = this->second(1);
  const Real & ot_min_y = other_box.first(1);
  const Real & ot_max_y = other_box.second(1);

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
  const Real & my_min_z = this->first(2);
  const Real & my_max_z = this->second(2);
  const Real & ot_min_z = other_box.first(2);
  const Real & ot_max_z = other_box.second(2);

  const bool z_int =
    is_between(my_min_z - abstol, ot_min_z, my_max_z + abstol) ||
    is_between(my_min_z - abstol, ot_max_z, my_max_z + abstol) ||
    is_between(ot_min_z - abstol, my_min_z, ot_max_z + abstol) ||
    is_between(ot_min_z - abstol, my_max_z, ot_max_z + abstol);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

bool BoundingBox::contains_point (const Point & p) const
{
  // Make local variables first to make things more clear in a moment
  Real my_min_x = this->first(0);
  Real my_max_x = this->second(0);
  bool x_int = is_between(my_min_x, p(0), my_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  Real my_min_y = this->first(1);
  Real my_max_y = this->second(1);
  bool y_int = is_between(my_min_y, p(1), my_max_y);

  intersection_true = intersection_true && y_int;
#endif


#if LIBMESH_DIM > 2
  Real my_min_z = this->first(2);
  Real my_max_z = this->second(2);
  bool z_int = is_between(my_min_z, p(2), my_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}


void BoundingBox::intersect_with (const BoundingBox & other_box)
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


void BoundingBox::union_with (const BoundingBox & other_box)
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



Real BoundingBox::signed_distance(const Point & p) const
{
  if (contains_point(p))
    {
      // Sign convention: if Point is inside the bbox, the distance is
      // negative. We then find the smallest distance to the different
      // sides of the box and return that.
      Real min_dist = std::numeric_limits<Real>::max();

      for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
        {
          min_dist = std::min(min_dist, std::abs(p(dir) - second(dir)));
          min_dist = std::min(min_dist, std::abs(p(dir) - first(dir)));
        }

      return -min_dist;
    }
  else // p is outside the box
    {
      Real dx[3] = {0., 0., 0.};

      // Compute distance "above"/"below" the box in each
      // direction. If the point is somewhere in between the (min,
      // max) values of the box, dx is 0.
      for (unsigned int dir=0; dir<LIBMESH_DIM; ++dir)
        {
          if (p(dir) > second(dir))
            dx[dir] = p(dir) - second(dir);
          else if (p(dir) < first(dir))
            dx[dir] = p(dir) - first(dir);
        }

      return std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }
}

} // namespace libMesh
