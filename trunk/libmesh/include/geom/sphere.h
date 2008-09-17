// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <cmath>


// Local includes
#include "surface.h"
#include "libmesh.h"


/**
 * This class defines a sphere.  It also computes coordinate
 * transformations between cartesian  \f$ (x, y, z) \f$
 * and spherical  \f$ (r, \theta, \phi) \f$ coordinates.
 * The spherical coordinates are valid in the ranges:  
 *
 * - \f$ 0 \le r      < \infty \f$
 * - \f$ 0 \le \theta < \pi \f$
 * - \f$ 0 \le \phi   < 2\pi \f$
 *
 * The coordinates are related as follows:
 * \f$ \phi \f$ is the angle in the xy plane
 * starting with 0. from the positive x axis,
 * \f$ \theta \f$ is measured against the positive
 * z axis.
   \verbatim

          \      | Z
           \theta|  
            \    |    .
             \   |   .
              \  |  .
               \ | .
                \|.
  ---------------+---------.---------
                /|\       .          Y
               /phi\     .  
              /  |  \   .  
             /   |   \ .  
            /.........\  
           /     |
        X /      
   \endverbatim
 *
 * \author Benjamin S. Kirk, Daniel Dreyer
 * \date 2002-2007
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
  Sphere (const Point& c, const Real r);

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
  void create_from_center_radius (const Point& c, const Real r);

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
  Real radius() const { return _rad; }

  /**
   * Returns the radius of the sphere as a writeable reference.
   */
  Real& radius() { return _rad; }

  /**
   * @returns the center of the sphere.
   */ 
  const Point& center() const { return _cent; }

  /**
   * @returns the center of the sphere.
   */ 
  Point& center() { return _cent; }

  /**
   * @returns the spherical coordinates for the
   * cartesian coordinates \p cart.
   */ 
  Point surface_coords (const Point& cart) const;

  /**
   * @returns the cartesian coordinates for the
   * spherical coordinates \p sph.
   */ 
  Point world_coords (const Point& sph) const;

  
private:

  
  /**
   * The center of the sphere.
   */
  Point _cent;

  /**
   * The radius of the sphere.
   */
  Real  _rad;
};



// ------------------------------------------------------------
// Sphere inline functions
inline
Point Sphere::surface_coords (const Point& cart) const
{
  // constant translation in the origin
  const Point c (cart-this->center());

  // phi: special care, so that it gives 0..2pi results
  const Real phi = std::atan2(c(1), c(0));

  return Point(/* radius */ c.size(),
	       /* theta  */ std::atan2( std::sqrt( c(0)*c(0) + c(1)*c(1) ), c(2) ),
	       /* phi    */ ( (phi < 0)  ?  2.*libMesh::pi+phi  :  phi ) );
}



inline
Point Sphere::world_coords (const Point& sph) const
{
  const Real r     = sph(0);
  const Real theta = sph(1);
  const Real phi   = sph(2);

  // constant translation out of the origin
  return Point (/* x */ r*std::sin(theta)*std::cos(phi) + this->center()(0),
		/* y */ r*std::sin(theta)*std::sin(phi) + this->center()(1),
		/* z */ r*std::cos(theta)               + this->center()(2));
}




#endif
