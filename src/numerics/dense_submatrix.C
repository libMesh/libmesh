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


// Local Includes
#include "libmesh/dense_submatrix.h"

namespace libMesh
{




// // ------------------------------------------------------------
// // Dense Matrix member functions
template<typename T>
void DenseSubMatrix<T>::left_multiply (const DenseMatrixBase<T> & M2)
{
  // (*this) <- M2 * M3
  // Where:
  // (*this) = (m x n),
  // M2      = (m x p),
  // M3      = (p x n)

  // M3 is a simply a copy of *this
  DenseSubMatrix<T> M3(*this);

  // Call the multiply function in the base class
  this->multiply(*this, M2, M3);
}



template<typename T>
void DenseSubMatrix<T>::right_multiply (const DenseMatrixBase<T> & M3)
{
  // (*this) <- M2 * M3
  // Where:
  // (*this) = (m x n),
  // M2      = (m x p),
  // M3      = (p x n)

  // M2 is simply a copy of *this
  DenseSubMatrix<T> M2(*this);

  // Call the multiply function in the base class
  this->multiply(*this, M2, M3);
}



//--------------------------------------------------------------
// Explicit instantiations
template class DenseSubMatrix<Real>;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class DenseSubMatrix<Complex>;
#endif

} // namespace libMesh
