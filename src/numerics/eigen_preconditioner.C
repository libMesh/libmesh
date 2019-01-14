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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// C++ includes

// Local Includes
#include "libmesh/eigen_preconditioner.h"
#include "libmesh/eigen_sparse_matrix.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/libmesh_common.h"

namespace libMesh
{

template <typename T>
void EigenPreconditioner<T>::apply(const NumericVector<T> & /* x */, NumericVector<T> & /* y */)
{
  libmesh_not_implemented();
}




template <typename T>
void EigenPreconditioner<T>::init ()
{
  libmesh_not_implemented();
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenPreconditioner<Number>;

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_EIGEN
