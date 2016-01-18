
// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PERIODIC_BOUNDARIES_H
#define LIBMESH_PERIODIC_BOUNDARIES_H

// ------------------------------------------------------------
// Periodic boundary conditions information

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes -----------------------------------
#include "libmesh/vector_value.h" // RealVectorValue

// C++ Includes   -----------------------------------
#include <map>

namespace libMesh {

// Forward Declarations
class Elem;
class PeriodicBoundaryBase;
class PointLocatorBase;

/**
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Note that std::map has no
 * virtual destructor, so downcasting here would be dangerous.
 */
class PeriodicBoundaries : public std::map<boundary_id_type, PeriodicBoundaryBase *>
{
public:
  PeriodicBoundaryBase * boundary(boundary_id_type id);

  const PeriodicBoundaryBase * boundary(boundary_id_type id) const;

  PeriodicBoundaries() {}

  ~PeriodicBoundaries();

  // The periodic neighbor of \p e in direction \p side, if it
  // exists.  NULL otherwise
  const Elem * neighbor(boundary_id_type boundary_id,
                        const PointLocatorBase & point_locator,
                        const Elem * e,
                        unsigned int side) const;
};

} // namespace libMesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARIES_H
