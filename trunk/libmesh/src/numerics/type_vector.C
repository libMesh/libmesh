// $Id: type_vector.C,v 1.7 2004-10-19 12:44:10 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include <iostream>
#include <iomanip> // for std::setw, std::setiosflags

// Local includes
#include "type_vector.h"




// ------------------------------------------------------------
// TypeVector<T> class member funcions
template <typename T>
TypeVector<T> TypeVector<T>::cross(const TypeVector<T>& p) const
{
  assert (DIM == 3);

  // |     i          j          k    |
  // |(*this)(0) (*this)(1) (*this)(2)|
  // |   p(0)       p(1)       p(2)   |
  
  return TypeVector<T>(  _coords[1]*p._coords[2] - _coords[2]*p._coords[1],
			-_coords[0]*p._coords[2] + _coords[2]*p._coords[0],
			 _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}



template <typename T>
TypeVector<T> TypeVector<T>::unit() const
{

  const Real length = size();
  
  assert (length != static_cast<Real>(0.));
  
#if DIM == 1
  return TypeVector<T>(_coords[0]/length);
#endif
  
#if DIM == 2 
  return TypeVector<T>(_coords[0]/length,
		       _coords[1]/length);
#endif
  
#if DIM == 3
  return TypeVector<T>(_coords[0]/length,
		       _coords[1]/length, 
		       _coords[2]/length);
#endif
  
}



template <typename T>
void TypeVector<T>::print() const
{
#if DIM == 1
  
  std::cout << "x=" << (*this)(0) << std::endl;
  
#endif
#if DIM == 2
  
  std::cout << "(x,y)=("
	    << std::setw(8) << (*this)(0) << ", "
	    << std::setw(8) << (*this)(1) << ")"
	    << std::endl;

#endif
#if DIM == 3
  
  std::cout <<  "(x,y,z)=("
	    << std::setw(8) << (*this)(0) << ", "
	    << std::setw(8) << (*this)(1) << ", "
	    << std::setw(8) << (*this)(2) << ")"
	    << std::endl;
#endif
}





template <typename T>
void TypeVector<T>::write_unformatted (std::ostream &out,
				       const bool newline) const
{
  assert (out);

  out << std::setiosflags(std::ios::showpoint)
      << (*this)(0) << " "
      << (*this)(1) << " "
      << (*this)(2) << " ";

  if (newline)
    out << '\n';      
}



template <>
bool TypeVector<Real>::operator < (const TypeVector<Real>& rhs) const
{
  if (*this == rhs)
    return false;


  
  if ((*this)(2) < rhs(2))
    return true;
    
  else
    if ((*this)(1) < rhs(1))
      return true;
  
      else
	if ((*this)(0) < rhs(0))
	  return true;
  
	else
	  return false;

  
  return false;
}



#ifdef USE_COMPLEX_NUMBERS
template <>
bool TypeVector<Complex>::operator < (const TypeVector<Complex>& rhs) const
{
  if (*this == rhs)
    return false;


  
  if ((*this)(2).real() < rhs(2).real() ||
      (*this)(2).imag() < rhs(2).imag())
    return true;
    
  else
    if ((*this)(1).real() < rhs(1).real() ||
	(*this)(1).imag() < rhs(1).imag())
      return true;
  
      else
	if ((*this)(0).real() < rhs(0).real() ||
	    (*this)(0).imag() < rhs(0).imag())
	  return true;
	
	else
	  return false;

  
  return false;
}



#endif



// ------------------------------------------------------------
// Explicit instantiations
template class TypeVector<Real>;

#ifdef USE_COMPLEX_NUMBERS
template class TypeVector<Complex>;
#endif
