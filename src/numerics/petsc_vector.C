// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/petsc_vector.h"

#ifdef LIBMESH_HAVE_PETSC

// libMesh includes
#include "libmesh/petsc_matrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/int_range.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/wrapped_petsc.h"

// TIMPI includes
#include "timpi/op_function.h"
#include "timpi/parallel_implementation.h"
#include "timpi/standard_type.h"

// C++ includes
#include <numeric> // std::iota

namespace libMesh
{
//-----------------------------------------------------------------------
// PetscVector members

template <typename T>
T PetscVector<T>::sum () const
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscScalar value=0.;

  ierr = VecSum (_vec, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<T>(value);
}


template <typename T>
Real PetscVector<T>::l1_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_1, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_2, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_INFINITY, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T> &
PetscVector<T>::operator += (const NumericVector<T> & v)
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert(this->closed());

  this->add(1., v);

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator -= (const NumericVector<T> & v)
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert(this->closed());

  this->add(-1., v);

  return *this;
}



template <typename T>
void PetscVector<T>::set (const numeric_index_type i, const T value)
{
  this->_restore_array();
  libmesh_assert_less (i, size());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscInt i_val = static_cast<PetscInt>(i);
  PetscScalar petsc_value = PS(value);

  std::scoped_lock lock(this->_numeric_vector_mutex);
  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::reciprocal()
{
  parallel_object_only();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // VecReciprocal has been in PETSc since at least 2.3.3 days
  ierr = VecReciprocal(_vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::conjugate()
{
  parallel_object_only();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // We just call the PETSc VecConjugate
  ierr = VecConjugate(_vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::add (const numeric_index_type i, const T value)
{
  this->_restore_array();
  libmesh_assert_less (i, size());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscInt i_val = static_cast<PetscInt>(i);
  PetscScalar petsc_value = PS(value);

  std::scoped_lock lock(this->_numeric_vector_mutex);
  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const T * v,
                                 const std::vector<numeric_index_type> & dof_indices)
{
  // If we aren't adding anything just return
  if (dof_indices.empty())
    return;

  this->_restore_array();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  const PetscInt * i_val = reinterpret_cast<const PetscInt *>(dof_indices.data());
  const PetscScalar * petsc_value = pPS(v);

  std::scoped_lock lock(this->_numeric_vector_mutex);
  ierr = VecSetValues (_vec, cast_int<PetscInt>(dof_indices.size()),
                       i_val, petsc_value, ADD_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T> & v_in,
                                 const SparseMatrix<T> & A_in)
{
  parallel_object_only();

  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // We shouldn't close() the matrix for you, as that would potentially modify the state of a const object.
  if (!A->closed())
    {
      libmesh_warning("Matrix A must be assembled before calling PetscVector::add_vector(v, A).\n"
                      "Please update your code, as this warning will become an error in a future release.");
      libmesh_deprecated();
      const_cast<PetscMatrix<T> *>(A)->close();
    }

  // The const_cast<> is not elegant, but it is required since PETSc
  // expects a non-const Mat.
  ierr = MatMultAdd(const_cast<PetscMatrix<T> *>(A)->mat(), v->_vec, _vec, _vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::add_vector_transpose (const NumericVector<T> & v_in,
                                           const SparseMatrix<T> & A_in)
{
  parallel_object_only();

  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // We shouldn't close() the matrix for you, as that would potentially modify the state of a const object.
  if (!A->closed())
    {
      libmesh_warning("Matrix A must be assembled before calling PetscVector::add_vector_transpose(v, A).\n"
                      "Please update your code, as this warning will become an error in a future release.");
      libmesh_deprecated();
      const_cast<PetscMatrix<T> *>(A)->close();
    }

  // The const_cast<> is not elegant, but it is required since PETSc
  // expects a non-const Mat.
  ierr = MatMultTransposeAdd(const_cast<PetscMatrix<T> *>(A)->mat(), v->_vec, _vec, _vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::add_vector_conjugate_transpose (const NumericVector<T> & v_in,
                                                     const SparseMatrix<T> & A_in)
{
  parallel_object_only();

  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  // We shouldn't close() the matrix for you, as that would potentially modify the state of a const object.
  if (!A->closed())
    {
      libmesh_warning("Matrix A must be assembled before calling PetscVector::add_vector_conjugate_transpose(v, A).\n"
                      "Please update your code, as this warning will become an error in a future release.");
      libmesh_deprecated();
      const_cast<PetscMatrix<T> *>(A)->close();
    }

  // Store a temporary copy since MatMultHermitianTransposeAdd doesn't seem to work
  // TODO: Find out why MatMultHermitianTransposeAdd doesn't work, might be a PETSc bug?
  std::unique_ptr<NumericVector<Number>> this_clone = this->clone();

  // The const_cast<> is not elegant, but it is required since PETSc
  // expects a non-const Mat.
  PetscErrorCode ierr = MatMultHermitianTranspose(const_cast<PetscMatrix<T> *>(A)->mat(), v->_vec, _vec);
  LIBMESH_CHKERR(ierr);

  // Add the temporary copy to the matvec result
  this->add(1., *this_clone);
}



template <typename T>
void PetscVector<T>::add (const T v_in)
{
  this->_get_array(false);

  for (numeric_index_type i=0; i<_local_size; i++)
    _values[i] += PetscScalar(v_in);
}



template <typename T>
void PetscVector<T>::add (const NumericVector<T> & v)
{
  parallel_object_only();

  this->add (1., v);
}



template <typename T>
void PetscVector<T>::add (const T a_in, const NumericVector<T> & v_in)
{
  parallel_object_only();

  this->_restore_array();

  // VecAXPY doesn't support &x==&y
  if (this == &v_in)
    {
      this->scale(a_in+1);
      return;
    }

  PetscScalar a = PS(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);
  v->_restore_array();

  libmesh_assert_equal_to (this->size(), v->size());

  PetscErrorCode ierr = VecAXPY(_vec, a, v->vec());
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);

  this->_is_closed = true;
}



template <typename T>
void PetscVector<T>::insert (const T * v,
                             const std::vector<numeric_index_type> & dof_indices)
{
  if (dof_indices.empty())
    return;

  this->_restore_array();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  PetscInt * idx_values = numeric_petsc_cast(dof_indices.data());
  std::scoped_lock lock(this->_numeric_vector_mutex);
  ierr = VecSetValues (_vec, cast_int<PetscInt>(dof_indices.size()),
                       idx_values, pPS(v), INSERT_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  parallel_object_only();

  this->_restore_array();

  PetscScalar factor = PS(factor_in);

  PetscErrorCode ierr = VecScale(_vec, factor);
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);
}

template <typename T>
NumericVector<T> & PetscVector<T>::operator *= (const NumericVector<T> & v)
{
  parallel_object_only();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  const PetscVector<T> * v_vec = cast_ptr<const PetscVector<T> *>(&v);

  ierr = VecPointwiseMult(_vec, _vec, v_vec->_vec);
  LIBMESH_CHKERR(ierr);

  return *this;
}

template <typename T>
NumericVector<T> & PetscVector<T>::operator /= (const NumericVector<T> & v)
{
  parallel_object_only();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  const PetscVector<T> * v_vec = cast_ptr<const PetscVector<T> *>(&v);

  ierr = VecPointwiseDivide(_vec, _vec, v_vec->_vec);
  LIBMESH_CHKERR(ierr);

  return *this;
}

template <typename T>
void PetscVector<T>::abs()
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr = VecAbs(_vec);
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);
}

template <typename T>
T PetscVector<T>::dot (const NumericVector<T> & v_in) const
{
  parallel_object_only();

  this->_restore_array();

  // Error flag
  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // Return value
  PetscScalar value=0.;

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<T>(value);
}

template <typename T>
T PetscVector<T>::indefinite_dot (const NumericVector<T> & v_in) const
{
  parallel_object_only();

  this->_restore_array();

  // Error flag
  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // Return value
  PetscScalar value=0.;

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecTDot(this->_vec, v->_vec, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<T>(value);
}


template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const T s_in)
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert(this->closed());

  PetscScalar s = PS(s_in);

  if (this->size() != 0)
    {
      PetscErrorCode ierr = VecSet(_vec, s);
      LIBMESH_CHKERR(ierr);

      libmesh_assert(this->comm().verify(int(this->type())));

      if (this->type() == GHOSTED)
        VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);
    }

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const NumericVector<T> & v_in)
{
  parallel_object_only();

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);

  *this = *v;

  return *this;
}



template <typename T>
PetscVector<T> &
PetscVector<T>::operator = (const PetscVector<T> & v)
{
  parallel_object_only();

  this->_restore_array();
  v._restore_array();

  libmesh_assert_equal_to (this->size(), v.size());
  libmesh_assert_equal_to (this->local_size(), v.local_size());
  libmesh_assert (v.closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  ierr = VecCopy (v._vec, this->_vec);
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);

  this->_is_closed = true;

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const std::vector<T> & v)
{
  parallel_object_only();

  this->_get_array(false);

  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    {
      numeric_index_type first = first_local_index();
      numeric_index_type last = last_local_index();
      for (numeric_index_type i=0; i<last-first; i++)
        _values[i] = PS(v[first + i]);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      for (numeric_index_type i=0; i<_local_size; i++)
        _values[i] = PS(v[i]);
    }

  // Make sure ghost dofs are up to date
  if (this->type() == GHOSTED)
    this->close();

  return *this;
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T> & v_local_in) const
{
  parallel_object_only();

  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * v_local = cast_ptr<PetscVector<T> *>(&v_local_in);

  libmesh_assert(v_local);
  // v_local_in should be closed
  libmesh_assert(v_local->closed());
  libmesh_assert_equal_to (v_local->size(), this->size());
  // 1) v_local_in is a large vector to hold the whole world
  // 2) v_local_in is a ghosted vector
  // 3) v_local_in is a parallel vector
  // Cases 2) and 3) should be scalable
  libmesh_assert(this->size()==v_local->local_size() || this->local_size()==v_local->local_size());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  if (v_local->type() == SERIAL && this->size() == v_local->local_size())
  {
    WrappedPetsc<VecScatter> scatter;
    ierr = VecScatterCreateToAll(_vec, scatter.get(), nullptr);
    LIBMESH_CHKERR(ierr);
    VecScatterBeginEnd(this->comm(), scatter, _vec, v_local->_vec, INSERT_VALUES, SCATTER_FORWARD);
  }
  // Two vectors have the same size, and we should just do a simple copy.
  // v_local could be either a parallel or ghosted vector
  else if (this->local_size() == v_local->local_size())
  {
    ierr = VecCopy(_vec,v_local->_vec);
    LIBMESH_CHKERR(ierr);
  }
  else
  {
    libmesh_error_msg("Vectors are inconsistent");
  }

  // Make sure ghost dofs are up to date
  // We do not call "close" here to save a global reduction
  if (v_local->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), v_local->_vec, INSERT_VALUES, SCATTER_FORWARD);
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T> & v_local_in,
                               const std::vector<numeric_index_type> & send_list) const
{
  parallel_object_only();

  libmesh_assert(this->comm().verify(int(this->type())));
  libmesh_assert(this->comm().verify(int(v_local_in.type())));

  // FIXME: Workaround for a strange bug at large-scale.
  // If we have ghosting, PETSc lets us just copy the solution, and
  // doing so avoids a segfault?
  if (v_local_in.type() == GHOSTED &&
      this->type() == PARALLEL)
    {
      v_local_in = *this;
      return;
    }

  // Normal code path begins here

  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * v_local = cast_ptr<PetscVector<T> *>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_equal_to (v_local->size(), this->size());
  libmesh_assert_less_equal (send_list.size(), v_local->size());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  const numeric_index_type n_sl =
    cast_int<numeric_index_type>(send_list.size());

  std::vector<PetscInt> idx(n_sl + this->local_size());
  for (numeric_index_type i=0; i<n_sl; i++)
    idx[i] = static_cast<PetscInt>(send_list[i]);
  for (auto i : make_range(this->local_size()))
    idx[n_sl+i] = i + this->first_local_index();

  // Create the index set & scatter objects
  WrappedPetsc<IS> is;
  PetscInt * idxptr = idx.empty() ? nullptr : idx.data();
  ierr = ISCreateGeneral(this->comm().get(), n_sl+this->local_size(),
                         idxptr, PETSC_USE_POINTER, is.get());
  LIBMESH_CHKERR(ierr);

  WrappedPetsc<VecScatter> scatter;
  ierr = VecScatterCreate(_vec,          is,
                          v_local->_vec, is,
                          scatter.get());
  LIBMESH_CHKERR(ierr);


  // Perform the scatter
  VecScatterBeginEnd(this->comm(), scatter, _vec, v_local->_vec, INSERT_VALUES, SCATTER_FORWARD);

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)
    v_local->close();
}



template <typename T>
void PetscVector<T>::localize (std::vector<T> & v_local,
                               const std::vector<numeric_index_type> & indices) const
{
  parallel_object_only();

  // Error code used to check the status of all PETSc function calls.
  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // Create a sequential destination Vec with the right number of entries on each proc.
  WrappedPetsc<Vec> dest;
  ierr = VecCreateSeq(PETSC_COMM_SELF, cast_int<PetscInt>(indices.size()), dest.get());
  LIBMESH_CHKERR(ierr);

  // Create an IS using the libmesh routine.  PETSc does not own the
  // IS memory in this case, it is automatically cleaned up by the
  // std::vector destructor.
  PetscInt * idxptr =
    indices.empty() ? nullptr : numeric_petsc_cast(indices.data());
  WrappedPetsc<IS> is;
  ierr = ISCreateGeneral(this->comm().get(), cast_int<PetscInt>(indices.size()), idxptr,
                         PETSC_USE_POINTER, is.get());
  LIBMESH_CHKERR(ierr);

  // Create the VecScatter object. "PETSC_NULL" means "use the identity IS".
  WrappedPetsc<VecScatter> scatter;
  ierr = VecScatterCreate(_vec,
                          /*src is=*/is,
                          /*dest vec=*/dest,
                          /*dest is=*/LIBMESH_PETSC_NULLPTR,
                          scatter.get());
  LIBMESH_CHKERR(ierr);

  // Do the scatter
  VecScatterBeginEnd(this->comm(), scatter, _vec, dest, INSERT_VALUES, SCATTER_FORWARD);

  // Get access to the values stored in dest.
  PetscScalar * values;
  ierr = VecGetArray (dest, &values);
  LIBMESH_CHKERR(ierr);

  // Store values into the provided v_local. Make sure there is enough
  // space reserved and then clear out any existing entries before
  // inserting.
  v_local.reserve(indices.size());
  v_local.clear();
  v_local.insert(v_local.begin(), values, values+indices.size());

  // We are done using it, so restore the array.
  ierr = VecRestoreArray (dest, &values);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::localize (const numeric_index_type first_local_idx,
                               const numeric_index_type last_local_idx,
                               const std::vector<numeric_index_type> & send_list)
{
  parallel_object_only();

  this->_restore_array();

  libmesh_assert_less_equal (send_list.size(), this->size());
  libmesh_assert_less_equal (last_local_idx+1, this->size());

  const numeric_index_type my_size       = this->size();
  const numeric_index_type my_local_size = (last_local_idx + 1 - first_local_idx);
  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // Don't bother for serial cases
  //  if ((first_local_idx == 0) &&
  //      (my_local_size == my_size))
  // But we do need to stay in sync for degenerate cases
  if (this->n_processors() == 1)
    return;


  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  PetscVector<T> parallel_vec(this->comm(), PARALLEL);

  parallel_vec.init (my_size, my_local_size, true, PARALLEL);


  // Copy part of *this into the parallel_vec
  {
    // Create idx, idx[i] = i+first_local_idx;
    std::vector<PetscInt> idx(my_local_size);
    std::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter objects
    WrappedPetsc<IS> is;
    ierr = ISCreateGeneral(this->comm().get(), my_local_size,
                           my_local_size ? idx.data() : nullptr, PETSC_USE_POINTER, is.get());
    LIBMESH_CHKERR(ierr);

    WrappedPetsc<VecScatter> scatter;
    ierr = VecScatterCreate(_vec,              is,
                            parallel_vec._vec, is,
                            scatter.get());
    LIBMESH_CHKERR(ierr);

    // Perform the scatter
    VecScatterBeginEnd(this->comm(), scatter, _vec, parallel_vec._vec, INSERT_VALUES, SCATTER_FORWARD);
  }

  // localize like normal
  parallel_vec.close();
  parallel_vec.localize (*this, send_list);
  this->close();
}



template <typename T>
void PetscVector<T>::localize (std::vector<T> & v_local) const
{
  parallel_object_only();

  this->_restore_array();

  // This function must be run on all processors at once
  parallel_object_only();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  const PetscInt n = this->size();
  const PetscInt nl = this->local_size();
  PetscScalar * values;

  v_local.clear();
  v_local.resize(n, 0.);

  ierr = VecGetArray (_vec, &values);
  LIBMESH_CHKERR(ierr);

  numeric_index_type ioff = first_local_index();

  for (PetscInt i=0; i<nl; i++)
    v_local[i+ioff] = static_cast<T>(values[i]);

  ierr = VecRestoreArray (_vec, &values);
  LIBMESH_CHKERR(ierr);

  this->comm().sum(v_local);
}



// Full specialization for Real datatypes
#ifdef LIBMESH_USE_REAL_NUMBERS

template <>
void PetscVector<Real>::localize_to_one (std::vector<Real> & v_local,
                                         const processor_id_type
                                         timpi_mpi_var(pid)) const
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  const PetscInt n  = size();
  PetscScalar * values;

  // only one processor
  if (n_processors() == 1)
    {
      v_local.resize(n);

      ierr = VecGetArray (_vec, &values);
      LIBMESH_CHKERR(ierr);

      for (PetscInt i=0; i<n; i++)
        v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
      LIBMESH_CHKERR(ierr);
    }

  // otherwise multiple processors
#ifdef LIBMESH_HAVE_MPI
  else
    {
      if (pid == 0) // optimized version for localizing to 0
        {
          WrappedPetsc<Vec> vout;
          WrappedPetsc<VecScatter> ctx;

          ierr = VecScatterCreateToZero(_vec, ctx.get(), vout.get());
          LIBMESH_CHKERR(ierr);

          VecScatterBeginEnd(this->comm(), ctx, _vec, vout, INSERT_VALUES, SCATTER_FORWARD);

          if (processor_id() == 0)
            {
              v_local.resize(n);

              ierr = VecGetArray (vout, &values);
              LIBMESH_CHKERR(ierr);

              for (PetscInt i=0; i<n; i++)
                v_local[i] = static_cast<Real>(values[i]);

              ierr = VecRestoreArray (vout, &values);
              LIBMESH_CHKERR(ierr);
            }
        }
      else
        {
          v_local.resize(n);

          numeric_index_type ioff = this->first_local_index();
          std::vector<Real> local_values (n, 0.);

          {
            ierr = VecGetArray (_vec, &values);
            LIBMESH_CHKERR(ierr);

            const PetscInt nl = local_size();
            for (PetscInt i=0; i<nl; i++)
              local_values[i+ioff] = static_cast<Real>(values[i]);

            ierr = VecRestoreArray (_vec, &values);
            LIBMESH_CHKERR(ierr);
          }


          MPI_Reduce (local_values.data(), v_local.data(), n,
                      Parallel::StandardType<Real>(),
                      Parallel::OpFunction<Real>::sum(), pid,
                      this->comm().get());
        }
    }
#endif // LIBMESH_HAVE_MPI
}

#endif


// Full specialization for Complex datatypes
#ifdef LIBMESH_USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize_to_one (std::vector<Complex> & v_local,
                                            const processor_id_type pid) const
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
  const PetscInt n  = size();
  const PetscInt nl = local_size();
  PetscScalar * values;


  v_local.resize(n);


  for (PetscInt i=0; i<n; i++)
    v_local[i] = 0.;

  // only one processor
  if (n_processors() == 1)
    {
      ierr = VecGetArray (_vec, &values);
      LIBMESH_CHKERR(ierr);

      for (PetscInt i=0; i<n; i++)
        v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
      LIBMESH_CHKERR(ierr);
    }

  // otherwise multiple processors
  else
    {
      numeric_index_type ioff = this->first_local_index();

      /* in here the local values are stored, acting as send buffer for MPI
       * initialize to zero, since we collect using MPI_SUM
       */
      std::vector<Real> real_local_values(n, 0.);
      std::vector<Real> imag_local_values(n, 0.);

      {
        ierr = VecGetArray (_vec, &values);
        LIBMESH_CHKERR(ierr);

        // provide my local share to the real and imag buffers
        for (PetscInt i=0; i<nl; i++)
          {
            real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
            imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
          }

        ierr = VecRestoreArray (_vec, &values);
        LIBMESH_CHKERR(ierr);
      }

      /* have buffers of the real and imaginary part of v_local.
       * Once MPI_Reduce() collected all the real and imaginary
       * parts in these std::vector<Real>, the values can be
       * copied to v_local
       */
      std::vector<Real> real_v_local(n);
      std::vector<Real> imag_v_local(n);

      // collect entries from other proc's in real_v_local, imag_v_local
      MPI_Reduce (real_local_values.data(),
                  real_v_local.data(), n,
                  Parallel::StandardType<Real>(),
                  Parallel::OpFunction<Real>::sum(),
                  pid, this->comm().get());

      MPI_Reduce (imag_local_values.data(),
                  imag_v_local.data(), n,
                  Parallel::StandardType<Real>(),
                  Parallel::OpFunction<Real>::sum(),
                  pid, this->comm().get());

      // copy real_v_local and imag_v_local to v_local
      for (PetscInt i=0; i<n; i++)
        v_local[i] = Complex(real_v_local[i], imag_v_local[i]);
    }
}

#endif



template <typename T>
void PetscVector<T>::pointwise_mult (const NumericVector<T> & vec1,
                                     const NumericVector<T> & vec2)
{
  parallel_object_only();

  this->_restore_array();

  // Convert arguments to PetscVector*.
  const PetscVector<T> * vec1_petsc = cast_ptr<const PetscVector<T> *>(&vec1);
  const PetscVector<T> * vec2_petsc = cast_ptr<const PetscVector<T> *>(&vec2);

  // Call PETSc function.
  PetscErrorCode ierr = VecPointwiseMult(_vec, vec1_petsc->vec(), vec2_petsc->vec());
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);

  this->_is_closed = true;
}

template <typename T>
void PetscVector<T>::pointwise_divide (const NumericVector<T> & vec1,
                                       const NumericVector<T> & vec2)
{
  parallel_object_only();

  this->_restore_array();

  // Convert arguments to PetscVector*.
  const PetscVector<T> * const vec1_petsc = cast_ptr<const PetscVector<T> *>(&vec1);
  const PetscVector<T> * const vec2_petsc = cast_ptr<const PetscVector<T> *>(&vec2);

  // Call PETSc function.
  PetscErrorCode ierr = VecPointwiseDivide(_vec, vec1_petsc->vec(), vec2_petsc->vec());
  LIBMESH_CHKERR(ierr);

  libmesh_assert(this->comm().verify(int(this->type())));

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);

  this->_is_closed = true;
}

template <typename T>
void PetscVector<T>::print_matlab (const std::string & name) const
{
  parallel_object_only();

  this->_restore_array();
  libmesh_assert (this->closed());

  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  WrappedPetsc<PetscViewer> petsc_viewer;
  ierr = PetscViewerCreate (this->comm().get(), petsc_viewer.get());
  LIBMESH_CHKERR(ierr);

  // Create an ASCII file containing the matrix
  // if a filename was provided.
  if (name != "")
    {
      ierr = PetscViewerASCIIOpen( this->comm().get(),
                                   name.c_str(),
                                   petsc_viewer.get());
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,7,0)
      ierr = PetscViewerSetFormat (petsc_viewer,
                                   PETSC_VIEWER_ASCII_MATLAB);
#else
      ierr = PetscViewerPushFormat (petsc_viewer,
                                    PETSC_VIEWER_ASCII_MATLAB);
#endif

      LIBMESH_CHKERR(ierr);

      ierr = VecView (_vec, petsc_viewer);
      LIBMESH_CHKERR(ierr);
    }

  // Otherwise the matrix will be dumped to the screen.
  else
    {

#if PETSC_VERSION_LESS_THAN(3,7,0)
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                   PETSC_VIEWER_ASCII_MATLAB);
#else
      ierr = PetscViewerPushFormat (PETSC_VIEWER_STDOUT_WORLD,
                                    PETSC_VIEWER_ASCII_MATLAB);
#endif

      LIBMESH_CHKERR(ierr);

      ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
      LIBMESH_CHKERR(ierr);
    }
}





template <typename T>
void PetscVector<T>::create_subvector(NumericVector<T> & subvector,
                                      const std::vector<numeric_index_type> & rows,
                                      const bool supplying_global_rows) const
{
  parallel_object_only();

  libmesh_error_msg_if(
      subvector.type() == GHOSTED,
      "We do not support scattering parallel information to ghosts for subvectors");

  this->_restore_array();

  // PETSc data structures
  PetscErrorCode ierr = static_cast<PetscErrorCode>(0);

  // Make sure the passed in subvector is really a PetscVector
  PetscVector<T> * petsc_subvector = cast_ptr<PetscVector<T> *>(&subvector);

  // If the petsc_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate PETSc routines to initialize it.
  if (!petsc_subvector->initialized())
    {
      libmesh_assert(petsc_subvector->_type == AUTOMATIC || petsc_subvector->_type == PARALLEL);

      if (supplying_global_rows)
        // Initialize the petsc_subvector to have enough space to hold
        // the entries which will be scattered into it.  Note: such an
        // init() function (where we let PETSc decide the number of local
        // entries) is not currently offered by the PetscVector
        // class.  Should we differentiate here between sequential and
        // parallel vector creation based on this->n_processors() ?
        ierr = VecCreateMPI(this->comm().get(),
                            PETSC_DECIDE,                    // n_local
                            cast_int<PetscInt>(rows.size()), // n_global
                            &(petsc_subvector->_vec));
      else
        ierr = VecCreateMPI(this->comm().get(),
                            cast_int<PetscInt>(rows.size()),
                            PETSC_DETERMINE,
                            &(petsc_subvector->_vec));
      LIBMESH_CHKERR(ierr);

      ierr = VecSetFromOptions (petsc_subvector->_vec);
      LIBMESH_CHKERR(ierr);

      // We created a parallel vector
      petsc_subvector->_type = PARALLEL;

      // Mark the subvector as initialized
      petsc_subvector->_is_initialized = true;
    }
  else
    {
      petsc_subvector->_restore_array();
    }

  std::vector<PetscInt> idx(rows.size());
  if (supplying_global_rows)
    std::iota (idx.begin(), idx.end(), 0);
  else
    {
      PetscInt start;
      ierr = VecGetOwnershipRange(petsc_subvector->_vec, &start, nullptr);
      LIBMESH_CHKERR(ierr);
      std::iota (idx.begin(), idx.end(), start);
    }

  // Construct index sets
  WrappedPetsc<IS> parent_is;
  ierr = ISCreateGeneral(this->comm().get(),
                         cast_int<PetscInt>(rows.size()),
                         numeric_petsc_cast(rows.data()),
                         PETSC_USE_POINTER,
                         parent_is.get());
  LIBMESH_CHKERR(ierr);

  WrappedPetsc<IS> subvector_is;
  ierr = ISCreateGeneral(this->comm().get(),
                         cast_int<PetscInt>(rows.size()),
                         idx.data(),
                         PETSC_USE_POINTER,
                         subvector_is.get());
  LIBMESH_CHKERR(ierr);

  // Construct the scatter object
  WrappedPetsc<VecScatter> scatter;
  ierr = VecScatterCreate(this->_vec,
                          parent_is,
                          petsc_subvector->_vec,
                          subvector_is,
                          scatter.get()); LIBMESH_CHKERR(ierr);

  // Actually perform the scatter
  VecScatterBeginEnd(this->comm(), scatter, this->_vec, petsc_subvector->_vec, INSERT_VALUES, SCATTER_FORWARD);

  petsc_subvector->_is_closed = true;
}



template <typename T>
void PetscVector<T>::_get_array(bool read_only) const
{
  libmesh_assert (this->initialized());

  bool initially_array_is_present = _array_is_present.load(std::memory_order_acquire);

  // If we already have a read/write array - and we're trying
  // to get a read only array - let's just use the read write
  if (initially_array_is_present && read_only && !_values_read_only)
    _read_only_values = _values;

  // If the values have already been retrieved and we're currently
  // trying to get a non-read only view (ie read/write) and the
  // values are currently read only... then we need to restore
  // the array first... and then retrieve it again.
  if (initially_array_is_present && !read_only && _values_read_only)
    {
      _restore_array();
      initially_array_is_present = false;
    }

  if (!initially_array_is_present)
    {
      std::scoped_lock lock(_petsc_get_restore_array_mutex);
      if (!_array_is_present.load(std::memory_order_relaxed))
        {
          PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
          if (this->type() != GHOSTED)
            {
              if (read_only)
                {
                  ierr = VecGetArrayRead(_vec, &_read_only_values);
                  _values_read_only = true;
                }
              else
                {
                  ierr = VecGetArray(_vec, &_values);
                  _values_read_only = false;
                }
              LIBMESH_CHKERR(ierr);
              _local_size = this->local_size();
            }
          else
            {
              ierr = VecGhostGetLocalForm (_vec,&_local_form);
              LIBMESH_CHKERR(ierr);

              if (read_only)
                {
                  ierr = VecGetArrayRead(_local_form, &_read_only_values);
                  _values_read_only = true;
                }
              else
                {
                  ierr = VecGetArray(_local_form, &_values);
                  _values_read_only = false;
                }
              LIBMESH_CHKERR(ierr);

              PetscInt my_local_size = 0;
              ierr = VecGetLocalSize(_local_form, &my_local_size);
              LIBMESH_CHKERR(ierr);
              _local_size = static_cast<numeric_index_type>(my_local_size);
            }

          { // cache ownership range
            PetscInt petsc_first=0, petsc_last=0;
            ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
            LIBMESH_CHKERR(ierr);
            _first = static_cast<numeric_index_type>(petsc_first);
            _last = static_cast<numeric_index_type>(petsc_last);
          }
          _array_is_present.store(true, std::memory_order_release);
        }
    }
}



template <typename T>
void PetscVector<T>::_restore_array() const
{
  libmesh_error_msg_if(_values_manually_retrieved,
                       "PetscVector values were manually retrieved but are being automatically restored!");

  libmesh_assert (this->initialized());
  if (_array_is_present.load(std::memory_order_acquire))
    {
      std::scoped_lock lock(_petsc_get_restore_array_mutex);
      if (_array_is_present.load(std::memory_order_relaxed))
        {
          PetscErrorCode ierr = static_cast<PetscErrorCode>(0);
          if (this->type() != GHOSTED)
            {
              if (_values_read_only)
                ierr = VecRestoreArrayRead (_vec, &_read_only_values);
              else
                ierr = VecRestoreArray (_vec, &_values);

              LIBMESH_CHKERR(ierr);
              _values = nullptr;
            }
          else
            {
              if (_values_read_only)
                ierr = VecRestoreArrayRead (_local_form, &_read_only_values);
              else
                ierr = VecRestoreArray (_local_form, &_values);

              LIBMESH_CHKERR(ierr);
              _values = nullptr;
              ierr = VecGhostRestoreLocalForm (_vec,&_local_form);
              LIBMESH_CHKERR(ierr);
              _local_form = nullptr;
              _local_size = 0;
            }
          _array_is_present.store(false, std::memory_order_release);
        }
    }
}

template <typename T>
std::unique_ptr<NumericVector<T>>
PetscVector<T>::get_subvector(const std::vector<numeric_index_type> & rows)
{
  // Construct index set
  WrappedPetsc<IS> parent_is;
  auto ierr = ISCreateGeneral(this->comm().get(),
                              cast_int<PetscInt>(rows.size()),
                              numeric_petsc_cast(rows.data()),
                              PETSC_USE_POINTER,
                              parent_is.get());
  LIBMESH_CHKERR(ierr);

  Vec subvec;
  ierr = VecGetSubVector(_vec, parent_is, &subvec);
  LIBMESH_CHKERR(ierr);

  return std::make_unique<PetscVector<T>>(subvec, this->comm());
}

template <typename T>
void
PetscVector<T>::restore_subvector(std::unique_ptr<NumericVector<T>> && subvector,
                                  const std::vector<numeric_index_type> & rows)
{
  auto * const petsc_subvector = cast_ptr<PetscVector<T> *>(subvector.get());

  // Construct index set
  WrappedPetsc<IS> parent_is;
  auto ierr = ISCreateGeneral(this->comm().get(),
                              cast_int<PetscInt>(rows.size()),
                              numeric_petsc_cast(rows.data()),
                              PETSC_USE_POINTER,
                              parent_is.get());
  LIBMESH_CHKERR(ierr);

  Vec subvec = petsc_subvector->vec();
  ierr = VecRestoreSubVector(_vec, parent_is, &subvec);
  LIBMESH_CHKERR(ierr);
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscVector<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
