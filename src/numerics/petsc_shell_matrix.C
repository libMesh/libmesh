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
#include "libmesh/petsc_shell_matrix.h"
#include "libmesh/petsc_matrix_shell_matrix.h"

namespace libMesh
{
template <typename T>
PetscShellMatrix<T>::~PetscShellMatrix()
{
  this->clear();
}

template <typename T>
void PetscShellMatrix<T>::vector_mult (NumericVector<T> & dest,
                                       const NumericVector<T> & arg) const
{
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);
  const PetscVector<T> & petsc_arg = cast_ref<const PetscVector<T> &>(arg);

  LibmeshPetscCall(MatMult(_mat, petsc_arg.vec(), petsc_dest.vec()));
}



template <typename T>
void PetscShellMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                           const NumericVector<T> & arg) const
{
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);
  const PetscVector<T> & petsc_arg = cast_ref<const PetscVector<T> &>(arg);

  LibmeshPetscCall(MatMultAdd(_mat, petsc_arg.vec(), petsc_dest.vec(), petsc_dest.vec()));
}


template <typename T>
void PetscShellMatrix<T>::clear ()
{
  if (this->initialized())
    {
      // If we encounter an error here, print a warning but otherwise
      // keep going since we may be recovering from an exception.
      PetscErrorCode ierr = MatDestroy (&_mat);
      if (ierr)
        libmesh_warning("Warning: MatDestroy returned a non-zero error code which we ignored.");

      this->_is_initialized = false;
    }
}

template <typename Obj>
void init_shell_mat(Obj & obj,
                    const numeric_index_type m,
                    const numeric_index_type n,
                    const numeric_index_type m_l,
                    const numeric_index_type n_l,
                    const numeric_index_type blocksize_in)
{
  // Clear initialized matrices
  if (obj.initialized())
    obj.clear();

  PetscInt m_global   = static_cast<PetscInt>(m);
  PetscInt n_global   = static_cast<PetscInt>(n);
  PetscInt m_local    = static_cast<PetscInt>(m_l);
  PetscInt n_local    = static_cast<PetscInt>(n_l);
  PetscInt blocksize  = static_cast<PetscInt>(blocksize_in);

  LibmeshPetscCall2(obj.comm(), MatCreate(obj.comm().get(), &obj._mat));
  LibmeshPetscCall2(obj.comm(), MatSetSizes(obj._mat, m_local, n_local, m_global, n_global));
  LibmeshPetscCall2(obj.comm(), MatSetBlockSize(obj._mat, blocksize));
  LibmeshPetscCall2(obj.comm(), MatSetType(obj._mat, MATSHELL));

  // Is prefix information available somewhere? Perhaps pass in the system name?
  LibmeshPetscCall2(obj.comm(), MatSetOptionsPrefix(obj._mat, ""));
  LibmeshPetscCall2(obj.comm(), MatSetFromOptions(obj._mat));
  LibmeshPetscCall2(obj.comm(), MatSetUp(obj._mat));
  LibmeshPetscCall2(obj.comm(), MatShellSetContext(obj._mat, &obj));

  obj._is_initialized = true;
}

template <typename Obj>
void init_shell_mat(Obj & obj)
{
  libmesh_assert(obj._dof_map);

  numeric_index_type my_m = obj._dof_map->n_dofs();
  numeric_index_type m_l = obj._dof_map->n_local_dofs();
  if (obj._omit_constrained_dofs)
    {
      my_m -= obj._dof_map->n_constrained_dofs();
      m_l -= obj._dof_map->n_local_constrained_dofs();
    }

  const numeric_index_type my_n = my_m;
  const numeric_index_type n_l  = m_l;

  init_shell_mat(obj, my_m, my_n, m_l, n_l, obj._dof_map->block_size());
}

template <typename T>
void PetscShellMatrix<T>::init ()
{
  init_shell_mat(*this);
}

template <typename T>
bool PetscShellMatrix<T>::initialized() const
{
  return _is_initialized;
}

template <typename T>
Mat PetscShellMatrix<T>::mat()
{
  libmesh_error_msg_if(!_mat, "A petsc shell matrix is not created yet. Please call init() first.");
  return _mat;
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscShellMatrix<Number>;
template void init_shell_mat(PetscShellMatrix<Number> & obj);
template void init_shell_mat(PetscMatrixShellMatrix<Number> & obj);
template void init_shell_mat(PetscShellMatrix<Number> & obj,
                             const numeric_index_type m,
                             const numeric_index_type n,
                             const numeric_index_type m_l,
                             const numeric_index_type n_l,
                             const numeric_index_type blocksize_in);
template void init_shell_mat(PetscMatrixShellMatrix<Number> & obj,
                             const numeric_index_type m,
                             const numeric_index_type n,
                             const numeric_index_type m_l,
                             const numeric_index_type n_l,
                             const numeric_index_type blocksize_in);

} // namespace libMesh

#endif
