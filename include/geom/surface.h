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



#ifndef LIBMESH_SURFACE_H
#define LIBMESH_SURFACE_H

// Local includes
#include "libmesh/point.h"

namespace libMesh
{

/**
 * The base class for all "surface" related geometric objects.  A
 * Surface is a two-dimensional object living in three-dimensional
 * space.  Examples of Surfaces are planes, hollow spheres, hollow
 * cylinders, etc...  This is a generic base class that describes the
 * useful functionality a surface will provide.  Specific derived
 * classes actually implement the functionality, so this class has
 * pure virtual members.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Base class for Plane and Sphere classes.
 */
class Surface
{
public:

  /**
   * Constructor.  Does nothing at the moment.
   */
  Surface () {}

  /**
   * Copy-constructor.
   */
  Surface (const Surface &) {}

  /**
   * Destructor.
   */
  virtual ~Surface () {}

  /**
   * \returns \p true if the point p is above the surface,
   * false otherwise.
   */
  virtual bool above_surface (const Point & p) const = 0;

  /**
   * \returns \p true if the point p is below the surface,
   * false otherwise.
   */
  virtual bool below_surface (const Point & p) const = 0;

  /**
   * \returns \p true if the point p is on the surface,
   * false otherwise.
   *
   * \note The definition of "on the surface" really means "very
   * close" to account for roundoff error.
   */
  virtual bool on_surface (const Point & p) const = 0;

  /**
   * \returns The closest point on the surface to point p.
   */
  virtual Point closest_point (const Point & p) const = 0;

  /**
   * \returns A unit vector normal to the surface at
   * point p.
   */
  virtual Point unit_normal (const Point & p) const = 0;

  /**
   * \returns The \p Point \p world_coords in the
   * surface's coordinate system.  \p world_coords
   * is in the world coordinate system.  This method
   * is not purely virtual, because there may be surfaces
   * that do not have their own coordinate system.  These
   * simply do not have to override this method.
   */
  virtual Point surface_coords (const Point & world_coords) const { return world_coords; }

  /**
   * \returns The world (cartesian) coordinates for the
   * surface coordinates \p surf_coords.  This method
   * is not purely virtual, because there may be surfaces
   * that do not have an own coordinate system.  These
   * simply do not have to override this method.
   */
  virtual Point world_coords (const Point & surf_coords) const { return surf_coords; }
};

} // namespace libMesh

#endif // LIBMESH_SURFACE_H
