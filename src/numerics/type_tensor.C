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

template <typename T>
bool TypeTensor<T>::is_hpd(Real rel_tol) const
{
  // Convenient reference to this object, to be used for matrix
  // operations.
  const auto & A = *this;

  // The norm of this matrix, needed for relative tolerance checks.
  auto A_norm = A.norm();

  // The zero matrix is only positive semi-definite
  if (A_norm == 0)
    return false;

  // An absolute tolerance (based on the user's provided relative
  // tolerance) to be used in the floating point comparisons below.
  auto abs_tol = rel_tol * A_norm;

  // Form the complex conjugate transpose of A
  TypeTensor<T> A_conjugate_transpose;
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      A_conjugate_transpose(i,j) = libmesh_conj(A(j,i));

  // Check if Hermitian
  if ((A - A_conjugate_transpose).norm() > abs_tol)
    return false;

  // If we made it here, then we are Hermitian, so now we just need to
  // check if we are positive-definite by checking that all principal
  // minors are positive. Note: the determinant of a Hermitian matrix
  // and all principal minors are real-valued. Since we already
  // checked that the matrix is Hermitian above, we don't bother to
  // check the complex parts are also zero now.
  if (std::real(A.det()) < -abs_tol)
    return false;

  // For 3x3, also check the upper 2x2 determinant
#if LIBMESH_DIM > 2
  if (std::real(A(0,0)*A(1,1) - A(0,1)*A(1,0)) < -abs_tol)
    return false;
#endif

  // For 3x3 and 2x2, also check the 1x1 determinant
#if LIBMESH_DIM > 1
  if (std::real(A(0,0)) < -abs_tol)
    return false;
#endif

  // If we made it here, then the matrix is Hermitian
  // positive-definite.
  return true;
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
template class LIBMESH_EXPORT TypeTensor<Real>;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class LIBMESH_EXPORT TypeTensor<Complex>;
#endif

} // namespace libMesh
