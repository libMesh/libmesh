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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/petsc_matrix_shell_matrix.h"

namespace libMesh
{

template <typename T>
void
PetscMatrixShellMatrix<T>::init(const numeric_index_type m,
                                const numeric_index_type n,
                                const numeric_index_type m_l,
                                const numeric_index_type n_l,
                                const numeric_index_type,
                                const numeric_index_type,
                                const numeric_index_type blocksize)
{
  init_shell_mat(*this, m, n, m_l, n_l, blocksize);

  this->set_context();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::init(ParallelType)
{
  init_shell_mat(*this);

  this->set_context();
}

template class LIBMESH_EXPORT PetscMatrixShellMatrix<Number>;

} // namespace libMesh

#endif
