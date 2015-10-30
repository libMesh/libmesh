// $Id: surface.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __surface_h__
#define __surface_h__

// C++ includes

// Local includes
#include "mesh_common.h"
#include "point.h"


/**
 * This class defines a surface.  A surface is a two-dimensional
 * object living in three-dimensional space.  Examples of surfaces
 * are planes, hollow spheres, hollow cylinders, etc...  This is 
 * a generic base class that describes the useful functionality
 * a surface will provide.  Specific derived classes actually implement
 * the functionality.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// Surface class definition
class Surface
{
public:
  
  /**
   * Constructor.  Does nothing at the moment.
   */
  Surface ();

  /**
   * Copy-constructor.
   */
  Surface (const Surface& other_surface);

  /**
   * Destructor.
   */
  virtual ~Surface ();

  /**
   * @returns true if the point p is above the surface,
   * false otherwise.
   */
  virtual bool above_surface (const Point& p) const;

  /**
   * @returns true if the point p is below the surface,
   * false otherwise.
   */
  virtual bool below_surface (const Point& p) const;

  /**
   * @returns true if the point p is on the surface,
   * false otherwise.  Note that the definition of on 
   * the surface really means "very close" to account 
   * for roundoff error.
   */
  virtual bool on_surface (const Point& p) const;

  /**
   * @returns the closest point on the surface to point p.
   */
  virtual Point closest_point (const Point& p) const;

  /**
   * @returns a unit vector normal to the surface at
   * point p.  
   */
  virtual Point unit_normal (const Point& p) const;

protected:

};

#endif
