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
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_matrix_shell_matrix.h"
#include "libmesh/petsc_matrix.h"

namespace libMesh
{


//-----------------------------------------------------------------------
// PetscMatrixBase members


// Constructor
template <typename T>
PetscMatrixBase<T>::PetscMatrixBase(const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _mat(nullptr),
  _destroy_mat_on_exit(true)
{}



// Constructor taking an existing Mat but not the responsibility
// for destroying it
template <typename T>
PetscMatrixBase<T>::PetscMatrixBase(Mat mat_in,
                                    const Parallel::Communicator & comm_in,
                                    const bool destroy_on_exit) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(destroy_on_exit)
{
  this->_mat = mat_in;
  this->_is_initialized = true;

  this->set_context();
}



// Destructor
template <typename T>
PetscMatrixBase<T>::~PetscMatrixBase()
{
  this->clear();
}



template <typename T>
void PetscMatrixBase<T>::clear () noexcept
{
  if ((this->initialized()) && (this->_destroy_mat_on_exit))
    {
      exceptionless_semiparallel_only();

      // If we encounter an error here, print a warning but otherwise
      // keep going since we may be recovering from an exception.
      PetscErrorCode ierr = MatDestroy (&_mat);
      if (ierr)
        libmesh_warning("Warning: MatDestroy returned a non-zero error code which we ignored.");

      this->_is_initialized = false;
    }
}

template <typename T>
void PetscMatrixBase<T>::set_destroy_mat_on_exit(bool destroy)
{
  this->_destroy_mat_on_exit = destroy;
}


template <typename T>
void PetscMatrixBase<T>::swap(PetscMatrixBase<T> & m_in)
{
  std::swap(_mat, m_in._mat);
  std::swap(_destroy_mat_on_exit, m_in._destroy_mat_on_exit);
}

template <typename T>
void PetscMatrixBase<T>::set_context()
{
  libmesh_assert(this->_mat);
  PetscContainer container;
  LibmeshPetscCall(PetscContainerCreate(this->comm().get(), &container));
  LibmeshPetscCall(PetscContainerSetPointer(container, this));
  LibmeshPetscCall(PetscObjectCompose((PetscObject)(Mat)this->_mat, "PetscMatrixCtx", (PetscObject)container));
  LibmeshPetscCall(PetscContainerDestroy(&container));
}

template <typename T>
PetscMatrixBase<T> * PetscMatrixBase<T>::get_context(Mat mat, const TIMPI::Communicator & comm)
{
  void * ctx;
  PetscContainer container;
  LibmeshPetscCall2(comm, PetscObjectQuery((PetscObject)mat, "PetscMatrixCtx", (PetscObject *)&container));
  if (!container)
    return nullptr;

  LibmeshPetscCall2(comm, PetscContainerGetPointer(container, &ctx));
  libmesh_assert(ctx);
  return static_cast<PetscMatrixBase<T> *>(ctx);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::m () const
{
  libmesh_assert (this->initialized());

  PetscInt petsc_m=0, petsc_n=0;

  LibmeshPetscCall(MatGetSize (this->_mat, &petsc_m, &petsc_n));

  return static_cast<numeric_index_type>(petsc_m);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::local_m () const
{
  libmesh_assert (this->initialized());

  PetscInt m = 0;

  LibmeshPetscCall(MatGetLocalSize (this->_mat, &m, NULL));

  return static_cast<numeric_index_type>(m);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::n () const
{
  libmesh_assert (this->initialized());

  PetscInt petsc_m=0, petsc_n=0;

  LibmeshPetscCall(MatGetSize (this->_mat, &petsc_m, &petsc_n));

  return static_cast<numeric_index_type>(petsc_n);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::local_n () const
{
  libmesh_assert (this->initialized());

  PetscInt n = 0;

  LibmeshPetscCall(MatGetLocalSize (this->_mat, NULL, &n));

  return static_cast<numeric_index_type>(n);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::row_start () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;

  LibmeshPetscCall(MatGetOwnershipRange(this->_mat, &start, &stop));

  return static_cast<numeric_index_type>(start);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::row_stop () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;

  LibmeshPetscCall(MatGetOwnershipRange(this->_mat, &start, &stop));

  return static_cast<numeric_index_type>(stop);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::col_start () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;

  LibmeshPetscCall(MatGetOwnershipRangeColumn(this->_mat, &start, &stop));

  return static_cast<numeric_index_type>(start);
}

template <typename T>
numeric_index_type PetscMatrixBase<T>::col_stop () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;

  LibmeshPetscCall(MatGetOwnershipRangeColumn(this->_mat, &start, &stop));

  return static_cast<numeric_index_type>(stop);
}

template <typename T>
void PetscMatrixBase<T>::close ()
{
  semiparallel_only();

  // BSK - 1/19/2004
  // strictly this check should be OK, but it seems to
  // fail on matrix-free matrices.  Do they falsely
  // state they are assembled?  Check with the developers...
  //   if (this->closed())
  //     return;

  MatAssemblyBeginEnd(this->comm(), this->_mat, MAT_FINAL_ASSEMBLY);
}

template <typename T>
bool PetscMatrixBase<T>::closed() const
{
  libmesh_assert (this->initialized());

  PetscBool assembled;

  LibmeshPetscCall(MatAssembled(this->_mat, &assembled));

  return (assembled == PETSC_TRUE);
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscMatrixBase<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
