// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs


// Local includes
#include "libmesh/tensor_value.h"
#include "libmesh/sphere.h"

namespace libMesh
{



// ------------------------------------------------------------
// Sphere class member functions
Sphere::Sphere () :
  _rad(-1.)
{
}



Sphere::Sphere (const Point & c,
                const Real r)
{
  libmesh_assert_greater (r, 0.);

  this->create_from_center_radius (c, r);
}



Sphere::Sphere (const Sphere & other_sphere) :
  Surface()
{
  this->create_from_center_radius (other_sphere.center(),
                                   other_sphere.radius());
}



Sphere::Sphere(const Point & pa,
               const Point & pb,
               const Point & pc,
               const Point & pd)
{
  Point pad = pa - pd;
  Point pbd = pb - pd;
  Point pcd = pc - pd;

  TensorValue<Real> T(pad,pbd,pcd);

  Real D = T.det();

  // The points had better not be coplanar
  libmesh_assert_greater (std::abs(D), 1e-12);

  Real e = 0.5*(pa.norm_sq() - pd.norm_sq());
  Real f = 0.5*(pb.norm_sq() - pd.norm_sq());
  Real g = 0.5*(pc.norm_sq() - pd.norm_sq());

  TensorValue<Real> T1(e,pad(1),pad(2),
                       f,pbd(1),pbd(2),
                       g,pcd(1),pcd(2));
  Real sx = T1.det()/D;

  TensorValue<Real> T2(pad(0),e,pad(2),
                       pbd(0),f,pbd(2),
                       pcd(0),g,pcd(2));
  Real sy = T2.det()/D;

  TensorValue<Real> T3(pad(0),pad(1),e,
                       pbd(0),pbd(1),f,
                       pcd(0),pcd(1),g);
  Real sz = T3.det()/D;

  Point c(sx,sy,sz);
  Real r = (c-pa).norm();

  this->create_from_center_radius(c,r);
}



Sphere::~Sphere ()
{
}



void Sphere::create_from_center_radius (const Point & c, const Real r)
{
  this->center() = c;
  this->radius() = r;

  libmesh_assert_greater (this->radius(), 0.);
}



bool Sphere::intersects (const Sphere & other_sphere) const
{
  return distance(other_sphere) < 0 ? true : false;
}



Real Sphere::distance (const Sphere & other_sphere) const
{
  libmesh_assert_greater ( this->radius(), 0. );
  libmesh_assert_greater ( other_sphere.radius(), 0. );

  const Real the_distance = (this->center() - other_sphere.center()).norm();

  return the_distance - (this->radius() + other_sphere.radius());
}



bool Sphere::above_surface (const Point & p) const
{
  libmesh_assert_greater (this->radius(), 0.);

  // create a vector from the center to the point.
  const Point w = p - this->center();

  if (w.norm() > this->radius())
    return true;

  return false;
}



bool Sphere::below_surface (const Point & p) const
{
  libmesh_assert_greater (this->radius(), 0.);

  return ( !this->above_surface (p) );
}



bool Sphere::on_surface (const Point & p) const
{
  libmesh_assert_greater (this->radius(), 0.);

  // Create a vector from the center to the point.
  const Point w = p - this->center();

  // if the size of that vector is the same as the radius() then
  // the point is on the surface.
  if (std::abs(w.norm() - this->radius()) < 1.e-10)
    return true;

  return false;
}



Point Sphere::closest_point (const Point & p) const
{
  libmesh_assert_greater (this->radius(), 0.);

  // get the normal from the surface in the direction
  // of p
  Point normal = this->unit_normal (p);

  // The closest point on the sphere is in the direction
  // of the normal a distance r from the center.
  const Point cp = this->center() + normal*this->radius();

  return cp;
}



Point Sphere::unit_normal (const Point & p) const
{
  libmesh_assert_greater (this->radius(), 0.);

  libmesh_assert_not_equal_to (p, this->center());

  // Create a vector from the center to the point
  Point n = p - this->center();

  return n.unit();
}

} // namespace libMesh
