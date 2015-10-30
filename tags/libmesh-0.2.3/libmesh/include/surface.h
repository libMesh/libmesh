// $Id: surface.h,v 1.5 2003-01-24 17:24:39 jwpeterson Exp $

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
#include "point.h"


/**
 * This class defines a surface.  A surface is a two-dimensional
 * object living in three-dimensional space.  Examples of surfaces
 * are planes, hollow spheres, hollow cylinders, etc...  This is 
 * a generic base class that describes the useful functionality
 * a surface will provide.  Specific derived classes actually implement
 * the functionality, so this class has pure virtual members.
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
  Surface () {};

  /**
   * Copy-constructor.
   */
  Surface (const Surface&) {};

  /**
   * Destructor.
   */
  virtual ~Surface () {};

  /**
   * @returns true if the point p is above the surface,
   * false otherwise.
   */
  virtual bool above_surface (const Point& p) const = 0;

  /**
   * @returns true if the point p is below the surface,
   * false otherwise.
   */
  virtual bool below_surface (const Point& p) const = 0;

  /**
   * @returns true if the point p is on the surface,
   * false otherwise.  Note that the definition of on 
   * the surface really means "very close" to account 
   * for roundoff error.
   */
  virtual bool on_surface (const Point& p) const = 0;

  /**
   * @returns the closest point on the surface to point p.
   */
  virtual Point closest_point (const Point& p) const = 0;

  /**
   * @returns a unit vector normal to the surface at
   * point p.  
   */
  virtual Point unit_normal (const Point& p) const = 0;

protected:

};

#endif
