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

// Local includes
#include "libmesh/face.h"
#include "libmesh/int_range.h"
#include "libmesh/point.h"

// C++ includes


namespace libMesh
{

Point Face::quasicircumcenter () const
{
  const Point p0 = this->point(0);

  // We have to use vertex-polygon area rather than real area if we
  // want this to be consistent for faces with curved edges
  Real total_area = 0;
  Point weighted_cc_sum;

  for (auto v : make_range(2u, this->n_vertices()))
    {
      const Point p1 = this->point(v-1);
      const Point p2 = this->point(v);

      const Point a = p1-p0;
      const Point b = p2-p0;
      const Point axb = a.cross(b);
      const Real norm2_axb = axb.norm_sq();
      const Real area = sqrt(norm2_axb)/2;

      const Point cc_minus_p0 = ((a.norm_sq() * b - b.norm_sq() * a).cross(axb)) / 2 / norm2_axb;

      weighted_cc_sum += area*cc_minus_p0;
      total_area += area;
    }

  weighted_cc_sum /= total_area;
  weighted_cc_sum += p0;

  return weighted_cc_sum;
}

} // namespace libMesh
