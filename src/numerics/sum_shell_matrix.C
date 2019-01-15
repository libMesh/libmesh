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
#include "libmesh/sum_shell_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/int_range.h"

namespace libMesh
{

template <typename T>
numeric_index_type SumShellMatrix<T>::m () const
{
  libmesh_assert(!matrices.empty());
  const numeric_index_type n_rows = matrices[0]->m();
#ifndef NDEBUG
  for (auto i : index_range(matrices))
    libmesh_assert_equal_to (matrices[i]->m(), n_rows);
#endif
  return n_rows;
}



template <typename T>
numeric_index_type SumShellMatrix<T>::n () const
{
  libmesh_assert(!matrices.empty());
  const numeric_index_type n_cols = matrices[0]->n();
#ifndef NDEBUG
  for (auto i : index_range(matrices))
    libmesh_assert_equal_to (matrices[i]->n(), n_cols);
#endif
  return n_cols;
}



template <typename T>
void SumShellMatrix<T>::vector_mult (NumericVector<T> & dest,
                                     const NumericVector<T> & arg) const
{
  dest.zero();
  this->vector_mult_add(dest,arg);
}



template <typename T>
void SumShellMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                         const NumericVector<T> & arg) const
{
  for (auto i : index_range(matrices))
    matrices[i]->vector_mult_add(dest, arg);
}



template <typename T>
void SumShellMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  std::unique_ptr<NumericVector<T>> a = dest.zero_clone();
  for (auto i : index_range(matrices))
    {
      matrices[i]->get_diagonal(*a);
      dest += *a;
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class SumShellMatrix<Number>;

} // namespace libMesh
