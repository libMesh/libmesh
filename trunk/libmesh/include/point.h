// $Id: point.h,v 1.9 2003-02-13 22:56:08 benkirk Exp $

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



#ifndef __point_h__
#define __point_h__

// C++ includes
#include <math.h>

// Local includes
#include "mesh_common.h"




/**
 * This class is based in no small part on the deal.II (http://www.dealii.org)
 * \p Point class.  This is the basic data type for spatial coordinates.
 * Using a \p Point allows for dimension-independed storage of things like
 * nodal values or solution gradients.  Accordingly, many of the operations
 * one would like to perform with spatial vectors are provided.
 */

class Point
{
 public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p DIM dimensions.
   */
  Point  (const Real x=0.,
	  const Real y=0.,
	  const Real z=0.);

  /**
   * Copy-constructor.
   */
  Point (const Point& p);

  /**
   * Destructor.
   */ 
  ~Point ();

  /**
   * Assign to a point without creating a temporary.
   */
  void assign (const Point &);

  /**
   * Return the \f$ i^{th} \f$ element of the point.
   */
  Real operator () (const unsigned int i) const;

  /**
   * Return a writeable reference to the \f$ i^{th} \f$ element of the point.
   */
  Real & operator () (const unsigned int i);
  
  /**
   * Add two points.
   */
  Point operator + (const Point &) const;

  /**
   * Add to this points.
   */
  const Point & operator += (const Point &);
  
  /**
   * Add to this point without creating a temporary.
   */
  void add (const Point &); 
  
  /**
   * Add a scaled value to this point without
   * creating a temporary.
   */
  void add_scaled (const Point &, const Real); 
  
  /**
   * Subtract two points.
   */
  Point operator - (const Point &) const;

  /**
   * Subtract from this point.
   */
  const Point & operator -= (const Point &);

  /**
   * Subtract from this point without creating a temporary.
   */
  void subtract (const Point &); 
  
  /**
   * Subtract a scaled value from this point without
   * creating a temporary.
   */
  void subtract_scaled (const Point &, const Real); 
  
  /**
   * Return the opposite of a point
   */
  Point operator - () const;
  
  /**
   * Multiply a point by a number, i.e. scale.
   */
  Point operator * (const Real) const;
  
  /**
   * Multiply this point by a number, i.e. scale.
   */
  const Point & operator *= (const Real);
  
  /**
   * Divide a point by a number, i.e. scale.
   */
  Point operator / (const Real) const;

  /**
   * Divide this point by a number, i.e. scale.
   */
  const Point & operator /= (const Real);

  /**
   * Multiply 2 points together, i.e. dot-product.
   */
  Real operator * (const Point &) const;

  /**
   * Multiply 2 points together, i.e. dot-product.
   */
  Point cross(const Point &) const;

  /**
   * Think of a point as a \p dim dimensional vector.  This
   * will return a unit vector aligned in that direction.
   */
  Point unit() const;

  /**
   * Reutrns the magnitude of the point, i.e. the square-root of the
   * sum of the elements squared.
   */
  Real size() const;

  /**
   * Zero the point in any dimension.
   */
  void zero();

  /**
   * @returns \p true if two points occupy the same
   * physical location in space.
   */
  bool operator == (const Point& rhs) const;
  
  /**
   * @returns \p true if this point is "less"
   * than another.  Useful for sorting.
   */
  bool operator < (const Point& rhs) const;

  /**
   * @returns a key associated with this point.  Useful for sorting.
   */
  unsigned int key() const;
  
  /**
   * Formatted print to \p std::cout.
   */
  void print() const;

  /**
   * Unformatted print to the stream \p out.  Simply prints the elements
   * of the point separated by spaces.
   */ 
  void write_unformatted (std::ostream &out) const;
    
 protected:

  /**
   * Reutrns the magnitude of the point squared, i.e. the square-root
   * of the sum of the elements squared.
   */
  Real size_sq() const;

  /**
   * The coordinates of the \p Point
   */
  Real _coords[DIM];

  /**
   * Make the derived class a friend
   */
  friend class Node;
};



//------------------------------------------------------
// Inline functions
inline
Point::Point (const Real x,
	      const Real y,
	      const Real z)
{
  _coords[0] = x;

  if (DIM > 1)
    {
      _coords[1] = y;

      if (DIM == 3)
	_coords[2] = z;
    }
}



inline
Point::Point (const Point& p)
{
  // copy the nodes from point p to me
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] = p._coords[i];
}



inline
Point::~Point ()
{
}



inline
void Point::assign (const Point &p)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] = p._coords[i];
}



inline
Real Point::operator () (const unsigned int i) const
{
  assert (i<3);

  if (i > (DIM-1))
    return 0.;
  
  return _coords[i];
}



inline
Real & Point::operator () (const unsigned int i)
{
  assert (i<DIM);
  
  return _coords[i];
}



inline
Point Point::operator + (const Point &p) const
{
 
#if DIM == 1
  return Point(_coords[0] + p._coords[0]);
#endif

#if DIM == 2 
  return Point(_coords[0] + p._coords[0],
	       _coords[1] + p._coords[1]);
#endif

#if DIM == 3
  return Point(_coords[0] + p._coords[0],
	       _coords[1] + p._coords[1],
	       _coords[2] + p._coords[2]);
#endif
	       
}



inline
const Point & Point::operator += (const Point &p)
{
  add (p);

  return *this;
}



inline
void Point::add (const Point &p)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] += p._coords[i];
}



inline
void Point::add_scaled (const Point &p, const Real factor)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] += factor*p._coords[i];
}



inline
Point Point::operator - (const Point &p) const
{

#if DIM == 1
  return Point(_coords[0] - p._coords[0]);
#endif

#if DIM == 2 
  return Point(_coords[0] - p._coords[0],
	       _coords[1] - p._coords[1]);
#endif

#if DIM == 3
  return Point(_coords[0] - p._coords[0],
	       _coords[1] - p._coords[1],
	       _coords[2] - p._coords[2]);
#endif

}



inline
const Point & Point::operator -= (const Point &p)
{
  subtract (p);

  return *this;
}



inline
void Point::subtract (const Point& p)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] -= p._coords[i];
}



inline
void Point::subtract_scaled (const Point &p, const Real factor)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] -= factor*p._coords[i];
}



inline
Point Point::operator - () const
{
  
#if DIM == 1
  return Point(-_coords[0]);
#endif

#if DIM == 2 
  return Point(-_coords[0],
	       -_coords[1]);
#endif

#if DIM == 3
  return Point(-_coords[0],
	       -_coords[1], 
	       -_coords[2]);
#endif
  
}



inline
Point Point::operator * (const Real factor) const
{

#if DIM == 1
  return Point(_coords[0]*factor);
#endif
  
#if DIM == 2 
  return Point(_coords[0]*factor,
	       _coords[1]*factor);
#endif
  
#if DIM == 3
  return Point(_coords[0]*factor,
	       _coords[1]*factor, 
	       _coords[2]*factor);
#endif
  
}




inline
const Point & Point::operator *= (const Real factor)
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] *= factor;

  return *this;
}




inline
Point Point::operator / (const Real factor) const
{
  assert (factor != static_cast<Real>(0.));

#if DIM == 1
  return Point(_coords[0]/factor);
#endif
  
#if DIM == 2 
  return Point(_coords[0]/factor,
	       _coords[1]/factor);
#endif
  
#if DIM == 3
  return Point(_coords[0]/factor,
	       _coords[1]/factor, 
	       _coords[2]/factor);
#endif
  
}




inline
const Point & Point::operator /= (const Real factor)
{
  assert (factor != static_cast<Real>(0.));

  for (unsigned int i=0; i<DIM; i++)
    _coords[i] /= factor;

  return *this;
}




inline
Real Point::operator * (const Point &p) const
{
#if DIM == 1
  return _coords[0]*p._coords[0];
#endif

#if DIM == 2
  return (_coords[0]*p._coords[0] +
	  _coords[1]*p._coords[1]);
#endif

#if DIM == 3
  return (_coords[0]*p._coords[0] +
	  _coords[1]*p._coords[1] +
	  _coords[2]*p._coords[2]);
#endif
}



inline
Real Point::size() const
{
  return sqrt(size_sq());  
}



inline
void Point::zero()
{
  for (unsigned int i=0; i<DIM; i++)
    _coords[i] = 0;
}



inline
Real Point::size_sq() const
{
  Real val = 0.;

  for (unsigned int i=0; i<DIM; i++)
    val += _coords[i]*_coords[i];

  return val;  
}



#endif
