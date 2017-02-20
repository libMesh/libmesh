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


// Local includes
#include "libmesh/cell.h"


namespace libMesh
{

BoundingBox Cell::loose_bounding_box () const
{
  // This might have curved sides, but it's definitely *not* curving
  // through 4-D space, so the full bounding box is just the merger of
  // the sides' bounding boxes.

  BoundingBox bbox = this->build_side_ptr(0)->loose_bounding_box();
  unsigned int my_n_sides = this->n_sides();
  for (unsigned s=1; s < my_n_sides; ++s)
    bbox.intersect(this->build_side_ptr(s)->loose_bounding_box());

  return bbox;
}

} // namespace libMesh
