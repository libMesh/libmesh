// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
PetscMatrixShellMatrix<T>::init(ParallelType libmesh_dbg_var(type))
{
#ifndef NDEBUG
  libmesh_assert(this->_dof_map);
  const auto m = this->_dof_map->n_dofs();
  const auto m_l = this->_dof_map->n_local_dofs();
  if (m != m_l)
    libmesh_assert(type != SERIAL);
#endif
  init_shell_mat(*this);
  this->set_context();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::zero()
{
  libmesh_error();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
PetscMatrixShellMatrix<T>::zero_clone() const
{
  libmesh_error();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
PetscMatrixShellMatrix<T>::clone() const
{
  libmesh_not_implemented();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::set(const numeric_index_type, const numeric_index_type, const T)
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::add(const numeric_index_type, const numeric_index_type, const T)
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::add_matrix(const DenseMatrix<T> &,
                                      const std::vector<numeric_index_type> &,
                                      const std::vector<numeric_index_type> &)
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::add_matrix(const DenseMatrix<T> &,
                                      const std::vector<numeric_index_type> &)
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::add(const T, const SparseMatrix<T> &)
{
  libmesh_error();
}

template <typename T>
T
PetscMatrixShellMatrix<T>::operator()(const numeric_index_type, const numeric_index_type) const
{
  libmesh_error();
}

template <typename T>
Real
PetscMatrixShellMatrix<T>::l1_norm() const
{
  libmesh_error();
}

template <typename T>
Real
PetscMatrixShellMatrix<T>::linfty_norm() const
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::print_personal(std::ostream &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::get_diagonal(NumericVector<T> &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::get_transpose(SparseMatrix<T> &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMatrixShellMatrix<T>::get_row(numeric_index_type,
                                   std::vector<numeric_index_type> &,
                                   std::vector<T> &) const
{
  libmesh_error();
}

template class LIBMESH_EXPORT PetscMatrixShellMatrix<Number>;

} // namespace libMesh

#endif
