// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/type_tensor.h"

namespace libMesh
{




// ------------------------------------------------------------
// TypeTensor<T> class member functions


template <typename T>
void TypeTensor<T>::print(std::ostream & os) const
{
#if LIBMESH_DIM == 1

  os << "x=" << (*this)(0,0) << std::endl;

#endif
#if LIBMESH_DIM == 2

  os << "(xx,xy)=("
     << std::setw(8) << (*this)(0,0) << ", "
     << std::setw(8) << (*this)(0,1) << ")"
     << std::endl;
  os << "(yx,yy)=("
     << std::setw(8) << (*this)(1,0) << ", "
     << std::setw(8) << (*this)(1,1) << ")"
     << std::endl;

#endif
#if LIBMESH_DIM == 3

  os <<  "(xx,xy,xz)=("
     << std::setw(8) << (*this)(0,0) << ", "
     << std::setw(8) << (*this)(0,1) << ", "
     << std::setw(8) << (*this)(0,2) << ")"
     << std::endl;
  os <<  "(yx,yy,yz)=("
     << std::setw(8) << (*this)(1,0) << ", "
     << std::setw(8) << (*this)(1,1) << ", "
     << std::setw(8) << (*this)(1,2) << ")"
     << std::endl;
  os <<  "(zx,zy,zz)=("
     << std::setw(8) << (*this)(2,0) << ", "
     << std::setw(8) << (*this)(2,1) << ", "
     << std::setw(8) << (*this)(2,2) << ")"
     << std::endl;
#endif
}





template <typename T>
void TypeTensor<T>::write_unformatted (std::ostream & out_stream,
                                       const bool newline) const
{
  libmesh_assert (out_stream);

  out_stream << std::setiosflags(std::ios::showpoint)
             << (*this)(0,0) << " "
             << (*this)(0,1) << " "
             << (*this)(0,2) << " ";
  if (newline)
    out_stream << '\n';

  out_stream << std::setiosflags(std::ios::showpoint)
             << (*this)(1,0) << " "
             << (*this)(1,1) << " "
             << (*this)(1,2) << " ";
  if (newline)
    out_stream << '\n';

  out_stream << std::setiosflags(std::ios::showpoint)
             << (*this)(2,0) << " "
             << (*this)(2,1) << " "
             << (*this)(2,2) << " ";
  if (newline)
    out_stream << '\n';
}



template <>
bool TypeTensor<Real>::operator < (const TypeTensor<Real> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      {
        if ((*this)(i,j) < rhs(i,j))
          return true;
        if ((*this)(i,j) > rhs(i,j))
          return false;
      }
  return false;
}



template <>
bool TypeTensor<Real>::operator > (const TypeTensor<Real> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      {
        if ((*this)(i,j) > rhs(i,j))
          return true;
        if ((*this)(i,j) < rhs(i,j))
          return false;
      }
  return false;
}



#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
bool TypeTensor<Complex>::operator < (const TypeTensor<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      {
        if ((*this)(i,j).real() < rhs(i,j).real())
          return true;
        if ((*this)(i,j).real() > rhs(i,j).real())
          return false;
        if ((*this)(i,j).imag() < rhs(i,j).imag())
          return true;
        if ((*this)(i,j).imag() > rhs(i,j).imag())
          return false;
      }
  return false;
}



template <>
bool TypeTensor<Complex>::operator > (const TypeTensor<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      {
        if ((*this)(i,j).real() > rhs(i,j).real())
          return true;
        if ((*this)(i,j).real() < rhs(i,j).real())
          return false;
        if ((*this)(i,j).imag() > rhs(i,j).imag())
          return true;
        if ((*this)(i,j).imag() < rhs(i,j).imag())
          return false;
      }
  return false;
}



#endif



// ------------------------------------------------------------
// Explicit instantiations
template class TypeTensor<Real>;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class TypeTensor<Complex>;
#endif

} // namespace libMesh
