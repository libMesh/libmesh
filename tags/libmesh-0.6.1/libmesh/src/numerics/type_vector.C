// $Id: type_vector.C,v 1.13 2007-10-21 20:48:52 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
void TypeVector<T>::print(std::ostream& os) const
{
#if DIM == 1
  
  os << "x=" << (*this)(0) << '\n';
  
#endif
#if DIM == 2
  
  os << "(x,y)=("
     << std::setw(8) << (*this)(0) << ", "
     << std::setw(8) << (*this)(1) << ")"
     << '\n';

#endif
#if DIM == 3
  
  os <<  "(x,y,z)=("
     << std::setw(8) << (*this)(0) << ", "
     << std::setw(8) << (*this)(1) << ", "
     << std::setw(8) << (*this)(2) << ")"
     << '\n';
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



template <typename T>
bool TypeVector<T>::operator < (const TypeVector<T>& rhs) const
{
  for (unsigned int i=0; i<DIM; i++)
    {
      if ((*this)(i) < rhs(i))
        return true;
      if ((*this)(i) > rhs(i))
        return false;
    }
  return false;
}



template <typename T>
bool TypeVector<T>::operator > (const TypeVector<T>& rhs) const
{
  for (unsigned int i=0; i<DIM; i++)
    {
      if ((*this)(i) > rhs(i))
        return true;
      if ((*this)(i) < rhs(i))
        return false;
    }
  return false;
}



#ifdef USE_COMPLEX_NUMBERS
template <>
bool TypeVector<Complex>::operator < (const TypeVector<Complex>& rhs) const
{
  for (unsigned int i=0; i<DIM; i++)
    {
      if ((*this)(i).real() < rhs(i).real())
        return true;
      if ((*this)(i).real() > rhs(i).real())
        return false;
      if ((*this)(i).imag() < rhs(i).imag())
        return true;
      if ((*this)(i).imag() > rhs(i).imag())
        return false;
    }
  return false;
}



template <>
bool TypeVector<Complex>::operator > (const TypeVector<Complex>& rhs) const
{
  for (unsigned int i=0; i<DIM; i++)
    {
      if ((*this)(i).real() > rhs(i).real())
        return true;
      if ((*this)(i).real() < rhs(i).real())
        return false;
      if ((*this)(i).imag() > rhs(i).imag())
        return true;
      if ((*this)(i).imag() < rhs(i).imag())
        return false;
    }
  return false;
}
#endif



// ------------------------------------------------------------
// Explicit instantiations
template class TypeVector<Real>;

#ifdef USE_COMPLEX_NUMBERS
template class TypeVector<Complex>;
#endif
