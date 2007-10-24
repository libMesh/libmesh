// $Id: plane.C,v 1.8 2003-09-02 18:02:42 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



// C++ includes

// Local includes
#include "plane.h"



// ------------------------------------------------------------
// Plane class member functions
Plane::Plane ()
{
}



Plane::Plane (const Point& p, 
	      const Point& n)
{
  create_from_point_normal (p, n);
}



Plane::Plane (const Point& p0, 
	      const Point& p1, 
	      const Point& p2)
{
  create_from_three_points (p0, p1, p2);
}



Plane::Plane (const Plane& other_plane) :
  Surface()
{
  create_from_point_normal(other_plane.point,
			   other_plane.normal);
}



Plane::~Plane ()
{
}



void Plane::create_from_point_normal (const Point& p, const Point& n)
{
  normal = n.unit();
  point  = p;
}



void Plane::create_from_three_points (const Point& p0,
				      const Point& p1,
				      const Point& p2)
{
  // Just use p0 for the point.
  point = p0;
  
  Point e0 = p1 - p0;
  Point e1 = p2 - p0;
  Point n  = e0.cross(e1);

  Point normal = n.unit(); 
}



void Plane::xy_plane (const Real zpos)
{
  const Point p (0., 0., zpos);
  const Point n (0., 0., 1.);

  point  = p;
  normal = n;
}



void Plane::xz_plane (const Real ypos)
{
  const Point p (0., ypos, 0.);
  const Point n (0., 1., 0.);

  point  = p;
  normal = n;
}



void Plane::yz_plane (const Real xpos)
{
  const Point p (xpos, 0., 0.);
  const Point n (1., 0., 0.);

  point  = p;
  normal = n;
}



bool Plane::above_surface (const Point& p) const 
{
  // Create a vector from the surface to point p;
  const Point w = p - point;

  // The point is above the surface if the projection
  // of that vector onto the normal is positive 
  const Real proj = w*normal;
  
  if (proj > 0.) 
    return true;

  return false;
}



bool Plane::below_surface (const Point& p) const 
{
  return ( !above_surface (p) );
}



bool Plane::on_surface (const Point& p) const 
{
  // Create a vector from the surface to point p;
  const Point w = p - point;

  // If the projection of that vector onto the
  // plane's normal is 0 then the point is in
  // the plane.
  if (w.size() < 1.e-10)
    return true;

  return false;
}



Point Plane::closest_point (const Point& p) const
{
  // Create a vector from the surface to point p;
  const Point w = p - point;

  // The closest point in the plane to point p
  // is in the negative direction normal direction
  // a distance w (dot) p.
  const Point cp = p - normal*(w*normal);
  
  return cp;
}



Point Plane::unit_normal (const Point&) const
{
  return normal;
}
