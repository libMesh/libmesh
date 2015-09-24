// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/periodic_boundaries.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/elem.h"
#include "libmesh/periodic_boundary.h"

namespace libMesh {

// ------------------------------------------------------------
// PeriodicBoundaries member functions


PeriodicBoundaries::~PeriodicBoundaries()
{
  for (std::map<boundary_id_type, PeriodicBoundaryBase*>::iterator it = begin(); it != end(); ++it)
    delete it->second;
}



PeriodicBoundaryBase* PeriodicBoundaries::boundary(boundary_id_type id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return NULL;
  return i->second;
}



const PeriodicBoundaryBase* PeriodicBoundaries::boundary(boundary_id_type id) const
{
  const_iterator i = this->find(id);
  if (i == this->end())
    return NULL;
  return i->second;
}




const Elem *PeriodicBoundaries::neighbor(boundary_id_type boundary_id,
                                         const PointLocatorBase& point_locator,
                                         const Elem* e,
                                         unsigned int side) const
{
  // Find a point on that side (and only that side)

  Point p = e->build_side(side)->centroid();

  const PeriodicBoundaryBase *b = this->boundary(boundary_id);
  libmesh_assert (b);
  p = b->get_corresponding_pos(p);

  return point_locator.operator()(p);
}

} // namespace libMesh





#endif // LIBMESH_ENABLE_PERIODIC
