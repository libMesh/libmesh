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



#ifndef LIBMESH_PLANE_H
#define LIBMESH_PLANE_H

// Local includes
#include "libmesh/surface.h"

// C++ includes

namespace libMesh
{

/**
 * This class defines a plane.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
class Plane : public Surface
{
public:

  /**
   * Dummy Constructor.
   */
  Plane ();

  /**
   * Constructs a plane containing point p with normal n.
   */
  Plane (const Point & p, const Point & n);

  /**
   * Constructs a plane containing the three points.  The
   * normal is determined in a counter-clockwise sense.  See
   * the create_from_three_points method for more details.
   */
  Plane (const Point & p0, const Point & p1, const Point & p2);

  /**
   * Copy-constructor.
   */
  Plane (const Plane & other_plane);

  /**
   * Destructor.  Does nothing at the moment.
   */
  ~Plane ();

  /**
   * Defines a plane containing point p with normal n.
   */
  void create_from_point_normal (const Point & p, const Point & n);

  /**
   * Defines a plane intersecting the three points
   * p0, p1, and p2.  The normal is constructed in a
   * counter-clockwise sense, i.e. (p1-p0)x(p2-p0);
   */
  void create_from_three_points (const Point & p0,
                                 const Point & p1,
                                 const Point & p2 );

  /**
   * Creates an XY plane located at z=zpos,
   */
  void xy_plane (const Real zpos=0.);

  /**
   * Creates an XZ plane located at y=ypos,
   */
  void xz_plane (const Real ypos=0.);

  /**
   * Creates an YZ plane located at x=xpos,
   */
  void yz_plane (const Real xpos=0.);

  /**
   * @returns true if the point p is above the surface,
   * false otherwise.
   */
  virtual bool above_surface (const Point & p) const libmesh_override;

  /**
   * @returns true if the point p is below the surface,
   * false otherwise.
   */
  virtual bool below_surface (const Point & p) const libmesh_override;

  /**
   * @returns true if the point p is on the surface,
   * false otherwise.  Note that the definition of on
   * the surface really means "very close" to account
   * for roundoff error.
   */
  virtual bool on_surface (const Point & p) const libmesh_override;

  /**
   * @returns the closest point on the surface to point p.
   */
  virtual Point closest_point (const Point & p) const libmesh_override;

  /**
   * @returns a unit vector normal to the surface at
   * point p.
   */
  virtual Point unit_normal (const Point & p) const libmesh_override;

  /**
   * @returns a point on the plane useful
   * for determining position
   */
  const Point & get_planar_point() const;


private:


  /**
   * Returns the normal for the plane.
   */
  const Point & normal () const;

  /**
   *  The plane is defined by a point and a normal.
   */
  Point _point;
  Point _normal;

};



// ------------------------------------------------------------
// Plane class inline members
inline const Point & Plane::normal () const
{
  return _normal;
}

} // namespace libMesh

#endif // LIBMESH_PLANE_H
