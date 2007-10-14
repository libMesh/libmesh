// $Id: sphere.C,v 1.14 2005-05-17 15:26:20 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <cmath> // for std::abs


// Local includes
#include "sphere.h"



// ------------------------------------------------------------
// Sphere class member functions
Sphere::Sphere () : 
  _rad(-1.)
{
}



Sphere::Sphere (const Point& c, 
		const Real   r)
{
  assert (r > 0.);

  this->create_from_center_radius (c, r);
}



Sphere::Sphere (const Sphere& other_sphere) :
  Surface()
{
  this->create_from_center_radius (other_sphere.center(),
				   other_sphere.radius());
}



Sphere::~Sphere ()
{
}



void Sphere::create_from_center_radius (const Point& c, const Real r)
{
  this->center() = c;
  this->radius() = r;

  assert (this->radius() > 0.);
}



bool Sphere::intersects (const Sphere& other_sphere) const
{
  assert ( this->radius() > 0. );
  assert ( other_sphere.radius() > 0. );

  const Real distance = (this->center() - other_sphere.center()).size();

  if (distance < (this->radius() + other_sphere.radius()) )
    return true;
  
  return false;
}



bool Sphere::above_surface (const Point& p) const 
{
  assert (this->radius() > 0.);

  // create a vector from the center to the point.
  const Point w = p - this->center();

  if (w.size() > this->radius())
    return true;

  return false;
}



bool Sphere::below_surface (const Point& p) const 
{
  assert (this->radius() > 0.);

  return ( !this->above_surface (p) );
}



bool Sphere::on_surface (const Point& p) const 
{
  assert (this->radius() > 0.);

  // Create a vector from the center to the point.
  const Point w = p - this->center();

  // if the size of that vector is the same as the radius() then
  // the point is on the surface.
  if (std::abs(w.size() - this->radius()) < 1.e-10)
    return true;

  return false;
}



Point Sphere::closest_point (const Point& p) const
{
  assert (this->radius() > 0.);

  // get the normal from the surface in the direction
  // of p
  Point normal = this->unit_normal (p);

  // The closest point on the sphere is in the direction
  // of the normal a distance r from the center.
  const Point cp = this->center() + normal*this->radius();
  
  return cp;
}



Point Sphere::unit_normal (const Point& p) const
{
  assert (this->radius() > 0.);

  assert ( !(p == this->center()) );
  
  // Create a vector from the center to the point
  Point n = p - this->center();

  return n.unit();
}
