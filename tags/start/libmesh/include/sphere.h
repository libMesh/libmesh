// $Id: sphere.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __sphere_h__
#define __sphere_h__

// C++ includes

// Local includes
#include "surface.h"


/**
 * This class defines a sphere.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// Sphere class definition
class Sphere : public Surface
{
public:
  
  /**
   * Dummy Constructor.
   */
  Sphere ();

  /**
   * Constructs a sphere of radius r centered at c.
   */
  Sphere (const Point& c, const real r);

  /**
   * Copy-constructor.
   */
  Sphere (const Sphere& other_sphere);

  /**
   * Destructor.  Does nothing at the moment.
   */
  ~Sphere ();

  /**
   * Defines a sphere of radius r centered at c.
   */
  void create_from_center_radius (const Point& c, const real r);

  /**
   * @returns true if other_sphere intersects this sphere,
   * false otherwise.
   */
  bool intersects (const Sphere& other_sphere) const;

  /**
   * @returns true if the point p is above the surface,
   * false otherwise.
   */
  bool above_surface (const Point& p) const;

  /**
   * @returns true if the point p is below the surface,
   * false otherwise.
   */
  bool below_surface (const Point& p) const;

  /**
   * @returns true if the point p is on the surface,
   * false otherwise.  Note that the definition of on 
   * the surface really means "very close" to account 
   * for roundoff error.
   */
  bool on_surface (const Point& p) const;

  /**
   * @returns the closest point on the surface to point p.
   */
  Point closest_point (const Point& p) const;

  /**
   * @returns a unit vector normal to the surface at
   * point p.  
   */
  Point unit_normal (const Point& p) const;

  /**
   * Returns the radius of the sphere.
   */
  real radius() const { return rad; };

  /**
   * @returns the center of the sphere.
   */ 
  const Point& center() const { return cent; };

  
private:

  
  /**
   * The center of the sphere.
   */
  Point cent;

  /**
   * The radius of the sphere.
   */
  real  rad;
};

#endif
