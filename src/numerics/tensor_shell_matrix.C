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



// Local includes
#include "libmesh/tensor_shell_matrix.h"

namespace libMesh
{

template <typename T>
void TensorShellMatrix<T>::vector_mult (NumericVector<T> & dest,
                                        const NumericVector<T> & arg) const
{
  dest = _v;
  dest.scale(_w.dot(arg));
}



template <typename T>
void TensorShellMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                            const NumericVector<T> & arg) const
{
  dest.add(_w.dot(arg),_v);
}



template <typename T>
void TensorShellMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  dest.pointwise_mult(_v,_w);
}



//------------------------------------------------------------------
// Explicit instantiations
template class TensorShellMatrix<Number>;

} // namespace libMesh
