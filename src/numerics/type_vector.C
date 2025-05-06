// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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




// Local includes
#include "libmesh/type_vector.h"

// C++ includes
#include <type_traits> // std::is_trivially_copyable


// Boost doesn't satisfy `is_trivially_copyable` for its float128
// wrapper, but surely it's memcpyable anyway?
#if !LIBMESH_DEFAULT_QUADRUPLE_PRECISION
static_assert(std::is_trivially_copyable<libMesh::TypeVector<libMesh::Real>>::value,
              "Someone made TypeVector non-TriviallyCopyable");
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
static_assert(std::is_trivially_copyable<libMesh::TypeVector<libMesh::Complex>>::value,
              "Someone made TypeVector non-TriviallyCopyable");
#endif
#endif // LIBMESH_DEFAULT_QUADRUPLE_PRECISION


namespace libMesh
{




// ------------------------------------------------------------
// TypeVector<T> class member functions


template <typename T>
void TypeVector<T>::write_unformatted (std::ostream & os,
                                       const bool newline) const
{
  os << std::setiosflags(std::ios::showpoint)
     << (*this)(0) << " "
     << (*this)(1) << " "
     << (*this)(2) << " ";

  if (newline)
    os << '\n';
}



template <typename T>
bool TypeVector<T>::operator < (const TypeVector<T> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      if ((*this)(i) < rhs(i))
        return true;
      if ((*this)(i) > rhs(i))
        return false;
    }
  return false;
}


template <typename T>
bool TypeVector<T>::operator <= (const TypeVector<T> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      if ((*this)(i) < rhs(i))
        return true;
      if ((*this)(i) > rhs(i))
        return false;
    }
  return true;
}



template <typename T>
bool TypeVector<T>::operator > (const TypeVector<T> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      if ((*this)(i) > rhs(i))
        return true;
      if ((*this)(i) < rhs(i))
        return false;
    }
  return false;
}


template <typename T>
bool TypeVector<T>::operator >= (const TypeVector<T> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    {
      if ((*this)(i) > rhs(i))
        return true;
      if ((*this)(i) < rhs(i))
        return false;
    }
  return true;
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
bool TypeVector<Complex>::operator < (const TypeVector<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
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
bool TypeVector<Complex>::operator <= (const TypeVector<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
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
  return true;
}



template <>
bool TypeVector<Complex>::operator > (const TypeVector<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
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



template <>
bool TypeVector<Complex>::operator >= (const TypeVector<Complex> & rhs) const
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
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
  return true;
}

#endif



template <>
auto
TypeVector<bool>::l1_norm() const
{
  bool ret{};
  for (const auto i : make_range(libmesh_dim))
    ret += _coords[i];

  return ret;
}



// ------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT TypeVector<Real>;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class LIBMESH_EXPORT TypeVector<Complex>;
#endif

} // namespace libMesh
