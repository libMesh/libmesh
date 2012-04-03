// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "sum_shell_matrix.h"
#include "numeric_vector.h"

namespace libMesh
{

template <typename T>
unsigned int SumShellMatrix<T>::m () const
{
  libmesh_assert(!matrices.empty());
  const unsigned int result = matrices[0]->m();
  for(unsigned int i=matrices.size(); i-->1; )
    {
      libmesh_assert(matrices[i]->m()==result);
    }
  return result;
}



template <typename T>
unsigned int SumShellMatrix<T>::n () const
{
  libmesh_assert(!matrices.empty());
  const unsigned int result = matrices[0]->n();
  for(unsigned int i=matrices.size(); i-->1; )
    {
      libmesh_assert(matrices[i]->n()==result);
    }
  return result;
}



template <typename T>
void SumShellMatrix<T>::vector_mult (NumericVector<T>& dest,
				     const NumericVector<T>& arg) const
{
  dest.zero();
  this->vector_mult_add(dest,arg);
}



template <typename T>
void SumShellMatrix<T>::vector_mult_add (NumericVector<T>& dest,
					 const NumericVector<T>& arg) const
{
  for(unsigned int i=matrices.size(); i-->0; )
    {
      matrices[i]->vector_mult_add(dest,arg);
    }
}



template <typename T>
void SumShellMatrix<T>::get_diagonal (NumericVector<T>& dest) const
{
  AutoPtr<NumericVector<T> > a = dest.clone();
  dest.zero();
  for(unsigned int i=matrices.size(); i-->0; )
    {
      matrices[i]->get_diagonal(*a);
      dest += *a;
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class SumShellMatrix<Number>;

} // namespace libMesh

