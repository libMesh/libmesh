// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

template <typename T>
void PetscShellMatrix<T>::vector_mult (NumericVector<T> & dest,
                                        const NumericVector<T> & arg) const
{
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);
  const PetscVector<T> & petsc_arg = cast_ref<const PetscVector<T> &>(arg);

  PetscErrorCode ierr = MatMult(_mat, petsc_arg.vec(), petsc_dest.vec());
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscShellMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                            const NumericVector<T> & arg) const
{
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);
  const PetscVector<T> & petsc_arg = cast_ref<const PetscVector<T> &>(arg);

  PetscErrorCode ierr = MatMultAdd(_mat, petsc_arg.vec(), petsc_dest.vec(), petsc_dest.vec());
  LIBMESH_CHKERR(ierr);
}


template <typename T>
void PetscShellMatrix<T>::clear ()
{
  if (this->initialized())
    {
      _mat.destroy();
      this->_is_initialized = false;
    }
}


template <typename T>
void PetscShellMatrix<T>::init ()
{
  libmesh_assert(this->_dof_map);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;


  const numeric_index_type my_m = this->_dof_map->n_dofs();
  const numeric_index_type my_n = my_m;
  const numeric_index_type n_l  = this->_dof_map->n_dofs_on_processor(this->processor_id());
  const numeric_index_type m_l  = n_l;


  PetscErrorCode ierr = 0;
  PetscInt m_global   = static_cast<PetscInt>(my_m);
  PetscInt n_global   = static_cast<PetscInt>(my_n);
  PetscInt m_local    = static_cast<PetscInt>(m_l);
  PetscInt n_local    = static_cast<PetscInt>(n_l);

  ierr = MatCreate(this->comm().get(), _mat.get());
  LIBMESH_CHKERR(ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
  LIBMESH_CHKERR(ierr);
  PetscInt blocksize  = static_cast<PetscInt>(this->_dof_map->block_size());
  ierr = MatSetBlockSize(_mat, blocksize);
  LIBMESH_CHKERR(ierr);

  ierr = MatSetType(_mat, MATSHELL);
  LIBMESH_CHKERR(ierr);

  // Is prefix information available somewhere? Perhaps pass in the system name?
  ierr = MatSetOptionsPrefix(_mat, "");
  LIBMESH_CHKERR(ierr);
  ierr = MatSetFromOptions(_mat);
  LIBMESH_CHKERR(ierr);
  ierr = MatSetUp(_mat);
  LIBMESH_CHKERR(ierr);
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
template class PetscShellMatrix<Number>;

} // namespace libMesh

#endif
