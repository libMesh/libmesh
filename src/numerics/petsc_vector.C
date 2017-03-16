// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local Includes
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/utility.h"

namespace libMesh
{

//-----------------------------------------------------------------------
// PetscVector members

template <typename T>
T PetscVector<T>::sum () const
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr=0;
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

  PetscErrorCode ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_1, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_2, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_INFINITY, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T> &
PetscVector<T>::operator += (const NumericVector<T> & v)
{
  this->_restore_array();
  libmesh_assert(this->closed());

  this->add(1., v);

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator -= (const NumericVector<T> & v)
{
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

  PetscErrorCode ierr=0;
  PetscInt i_val = static_cast<PetscInt>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::reciprocal()
{
  PetscErrorCode ierr = 0;

  // VecReciprocal has been in PETSc since at least 2.3.3 days
  ierr = VecReciprocal(_vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::conjugate()
{
  PetscErrorCode ierr = 0;

  // We just call the PETSc VecConjugate
  ierr = VecConjugate(_vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::add (const numeric_index_type i, const T value)
{
  this->_restore_array();
  libmesh_assert_less (i, size());

  PetscErrorCode ierr=0;
  PetscInt i_val = static_cast<PetscInt>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const T * v,
                                 const std::vector<numeric_index_type> & dof_indices)
{
  // If we aren't adding anything just return
  if(dof_indices.empty())
    return;

  this->_restore_array();

  PetscErrorCode ierr=0;
  const PetscInt * i_val = reinterpret_cast<const PetscInt *>(&dof_indices[0]);
  const PetscScalar * petsc_value = static_cast<const PetscScalar *>(v);

  ierr = VecSetValues (_vec, cast_int<PetscInt>(dof_indices.size()),
                       i_val, petsc_value, ADD_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T> & V_in,
                                 const SparseMatrix<T> & A_in)
{
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * V = cast_ptr<const PetscVector<T> *>(&V_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  PetscErrorCode ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultAdd(const_cast<PetscMatrix<T> *>(A)->mat(), V->_vec, _vec, _vec);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::add_vector_transpose (const NumericVector<T> & V_in,
                                           const SparseMatrix<T> & A_in)
{
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * V = cast_ptr<const PetscVector<T> *>(&V_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  PetscErrorCode ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultTransposeAdd(const_cast<PetscMatrix<T> *>(A)->mat(), V->_vec, _vec, _vec);
  LIBMESH_CHKERR(ierr);
}


#if PETSC_VERSION_LESS_THAN(3,1,0)
template <typename T>
void PetscVector<T>::add_vector_conjugate_transpose (const NumericVector<T> &,
                                                     const SparseMatrix<T> &)
{
  libmesh_error_msg("MatMultHermitianTranspose was introduced in PETSc 3.1.0," \
                    << "No one has made it backwards compatible with older " \
                    << "versions of PETSc so far.");
}

#else

template <typename T>
void PetscVector<T>::add_vector_conjugate_transpose (const NumericVector<T> & V_in,
                                                     const SparseMatrix<T> & A_in)
{
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T> * V = cast_ptr<const PetscVector<T> *>(&V_in);
  const PetscMatrix<T> * A = cast_ptr<const PetscMatrix<T> *>(&A_in);

  A->close();

  // Store a temporary copy since MatMultHermitianTransposeAdd doesn't seem to work
  // TODO: Find out why MatMultHermitianTransposeAdd doesn't work, might be a PETSc bug?
  UniquePtr< NumericVector<Number> > this_clone = this->clone();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  PetscErrorCode ierr = MatMultHermitianTranspose(const_cast<PetscMatrix<T> *>(A)->mat(), V->_vec, _vec);
  LIBMESH_CHKERR(ierr);

  // Add the temporary copy to the matvec result
  this->add(1., *this_clone);
}
#endif


template <typename T>
void PetscVector<T>::add (const T v_in)
{
  this->_get_array(false);

  for (numeric_index_type i=0; i<_local_size; i++)
    _values[i] += v_in;
}



template <typename T>
void PetscVector<T>::add (const NumericVector<T> & v)
{
  this->add (1., v);
}



template <typename T>
void PetscVector<T>::add (const T a_in, const NumericVector<T> & v_in)
{
  this->_restore_array();

  PetscErrorCode ierr = 0;
  PetscScalar a = static_cast<PetscScalar>(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);
  v->_restore_array();

  libmesh_assert_equal_to (this->size(), v->size());

  if(this->type() != GHOSTED)
    {
      ierr = VecAXPY(_vec, a, v->_vec);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      Vec loc_vec;
      Vec v_loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostGetLocalForm (v->_vec,&v_loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecAXPY(loc_vec, a, v_loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecGhostRestoreLocalForm (v->_vec,&v_loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
    }
}



template <typename T>
void PetscVector<T>::insert (const T * v,
                             const std::vector<numeric_index_type> & dof_indices)
{
  if (dof_indices.empty())
    return;

  this->_restore_array();

  PetscErrorCode ierr=0;
  PetscInt * idx_values = numeric_petsc_cast(&dof_indices[0]);
  ierr = VecSetValues (_vec, dof_indices.size(), idx_values, v, INSERT_VALUES);
  LIBMESH_CHKERR(ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  this->_restore_array();

  PetscErrorCode ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);

  if(this->type() != GHOSTED)
    {
      ierr = VecScale(_vec, factor);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecScale(loc_vec, factor);
      LIBMESH_CHKERR(ierr);

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
    }
}

template <typename T>
NumericVector<T> & PetscVector<T>::operator /= (NumericVector<T> & v)
{
  PetscErrorCode ierr = 0;

  const PetscVector<T> * v_vec = cast_ptr<const PetscVector<T> *>(&v);

  ierr = VecPointwiseDivide(_vec, _vec, v_vec->_vec);
  LIBMESH_CHKERR(ierr);

  return *this;
}

template <typename T>
void PetscVector<T>::abs()
{
  this->_restore_array();

  PetscErrorCode ierr = 0;

  if(this->type() != GHOSTED)
    {
      ierr = VecAbs(_vec);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecAbs(loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
    }
}

template <typename T>
T PetscVector<T>::dot (const NumericVector<T> & V) const
{
  this->_restore_array();

  // Error flag
  PetscErrorCode ierr = 0;

  // Return value
  PetscScalar value=0.;

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&V);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<T>(value);
}

template <typename T>
T PetscVector<T>::indefinite_dot (const NumericVector<T> & V) const
{
  this->_restore_array();

  // Error flag
  PetscErrorCode ierr = 0;

  // Return value
  PetscScalar value=0.;

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&V);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecTDot(this->_vec, v->_vec, &value);
  LIBMESH_CHKERR(ierr);

  return static_cast<T>(value);
}


template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const T s_in)
{
  this->_restore_array();
  libmesh_assert(this->closed());

  PetscErrorCode ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)
    {
      if(this->type() != GHOSTED)
        {
          ierr = VecSet(_vec, s);
          LIBMESH_CHKERR(ierr);
        }
      else
        {
          Vec loc_vec;
          ierr = VecGhostGetLocalForm (_vec,&loc_vec);
          LIBMESH_CHKERR(ierr);

          ierr = VecSet(loc_vec, s);
          LIBMESH_CHKERR(ierr);

          ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
          LIBMESH_CHKERR(ierr);
        }
    }

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const NumericVector<T> & v_in)
{
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T> * v = cast_ptr<const PetscVector<T> *>(&v_in);

  *this = *v;

  return *this;
}



template <typename T>
PetscVector<T> &
PetscVector<T>::operator = (const PetscVector<T> & v)
{
  this->_restore_array();
  v._restore_array();

  libmesh_assert_equal_to (this->size(), v.size());
  libmesh_assert_equal_to (this->local_size(), v.local_size());
  libmesh_assert (v.closed());

  PetscErrorCode ierr = 0;

  if (((this->type()==PARALLEL) && (v.type()==GHOSTED)) ||
      ((this->type()==GHOSTED) && (v.type()==PARALLEL)) ||
      ((this->type()==GHOSTED) && (v.type()==SERIAL))   ||
      ((this->type()==SERIAL) && (v.type()==GHOSTED)))
    {
      /* Allow assignment of a ghosted to a parallel vector since this
         causes no difficulty.  See discussion in libmesh-devel of
         June 24, 2010.  */
      ierr = VecCopy (v._vec, this->_vec);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      /* In all other cases, we assert that both vectors are of equal
         type.  */
      libmesh_assert_equal_to (this->_type, v._type);

      if (v.size() != 0)
        {
          if(this->type() != GHOSTED)
            {
              ierr = VecCopy (v._vec, this->_vec);
              LIBMESH_CHKERR(ierr);
            }
          else
            {
              Vec loc_vec;
              Vec v_loc_vec;
              ierr = VecGhostGetLocalForm (_vec,&loc_vec);
              LIBMESH_CHKERR(ierr);
              ierr = VecGhostGetLocalForm (v._vec,&v_loc_vec);
              LIBMESH_CHKERR(ierr);

              ierr = VecCopy (v_loc_vec, loc_vec);
              LIBMESH_CHKERR(ierr);

              ierr = VecGhostRestoreLocalForm (v._vec,&v_loc_vec);
              LIBMESH_CHKERR(ierr);
              ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
              LIBMESH_CHKERR(ierr);
            }
        }
    }

  close();

  return *this;
}



template <typename T>
NumericVector<T> &
PetscVector<T>::operator = (const std::vector<T> & v)
{
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
        _values[i] =  static_cast<PetscScalar>(v[first + i]);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      for (numeric_index_type i=0; i<_local_size; i++)
        _values[i] = static_cast<PetscScalar>(v[i]);
    }

  // Make sure ghost dofs are up to date
  if (this->type() == GHOSTED)
    this->close();

  return *this;
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T> & v_local_in) const
{
  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * v_local = cast_ptr<PetscVector<T> *>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_equal_to (v_local->size(), this->size());

  PetscErrorCode ierr = 0;
  const PetscInt n = this->size();

  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<PetscInt> idx(n); Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateLibMesh(this->comm().get(), n, &idx[0], PETSC_USE_POINTER, &is);
  LIBMESH_CHKERR(ierr);

  ierr = VecScatterCreate(_vec,          is,
                          v_local->_vec, is,
                          &scatter);
  LIBMESH_CHKERR(ierr);

  // Perform the scatter
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

  // Clean up
  ierr = LibMeshISDestroy (&is);
  LIBMESH_CHKERR(ierr);

  ierr = LibMeshVecScatterDestroy(&scatter);
  LIBMESH_CHKERR(ierr);

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)
    v_local->close();
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T> & v_local_in,
                               const std::vector<numeric_index_type> & send_list) const
{
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

  PetscErrorCode ierr=0;
  const numeric_index_type n_sl =
    cast_int<numeric_index_type>(send_list.size());

  IS is;
  VecScatter scatter;

  std::vector<PetscInt> idx(n_sl + this->local_size());

  for (numeric_index_type i=0; i<n_sl; i++)
    idx[i] = static_cast<PetscInt>(send_list[i]);
  for (numeric_index_type i = 0; i != this->local_size(); ++i)
    idx[n_sl+i] = i + this->first_local_index();

  // Create the index set & scatter object
  if (idx.empty())
    ierr = ISCreateLibMesh(this->comm().get(),
                           n_sl+this->local_size(), PETSC_NULL, PETSC_USE_POINTER, &is);
  else
    ierr = ISCreateLibMesh(this->comm().get(),
                           n_sl+this->local_size(), &idx[0], PETSC_USE_POINTER, &is);
  LIBMESH_CHKERR(ierr);

  ierr = VecScatterCreate(_vec,          is,
                          v_local->_vec, is,
                          &scatter);
  LIBMESH_CHKERR(ierr);


  // Perform the scatter
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

  // Clean up
  ierr = LibMeshISDestroy (&is);
  LIBMESH_CHKERR(ierr);

  ierr = LibMeshVecScatterDestroy(&scatter);
  LIBMESH_CHKERR(ierr);

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)
    v_local->close();
}



template <typename T>
void PetscVector<T>::localize (std::vector<T> & v_local,
                               const std::vector<numeric_index_type> & indices) const
{
  // Error code used to check the status of all PETSc function calls.
  PetscErrorCode ierr = 0;

  // Create a sequential destination Vec with the right number of entries on each proc.
  Vec dest;
  ierr = VecCreateSeq(PETSC_COMM_SELF, indices.size(), &dest);
  LIBMESH_CHKERR(ierr);

  // Create an IS using the libmesh routine.  PETSc does not own the
  // IS memory in this case, it is automatically cleaned up by the
  // std::vector destructor.
  IS is;
  ierr = ISCreateLibMesh(this->comm().get(),
                         indices.size(),
                         numeric_petsc_cast(&indices[0]),
                         PETSC_USE_POINTER,
                         &is);
  LIBMESH_CHKERR(ierr);

  // Create the VecScatter object.  "NULL" means "use the identity IS".
  VecScatter scatter;
  ierr = VecScatterCreate(_vec,
                          /*src is=*/is,
                          /*dest vec=*/dest,
                          /*dest is=*/NULL,
                          &scatter);
  LIBMESH_CHKERR(ierr);

  // Do the scatter
  ierr = VecScatterBegin(scatter, _vec, dest, INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

  ierr = VecScatterEnd(scatter, _vec, dest, INSERT_VALUES, SCATTER_FORWARD);
  LIBMESH_CHKERR(ierr);

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

  // Clean up PETSc data structures.
  ierr = LibMeshVecScatterDestroy(&scatter);
  LIBMESH_CHKERR(ierr);

  ierr = LibMeshISDestroy (&is);
  LIBMESH_CHKERR(ierr);

  ierr = LibMeshVecDestroy(&dest);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::localize (const numeric_index_type first_local_idx,
                               const numeric_index_type last_local_idx,
                               const std::vector<numeric_index_type> & send_list)
{
  this->_restore_array();

  libmesh_assert_less_equal (send_list.size(), this->size());
  libmesh_assert_less_equal (last_local_idx+1, this->size());

  const numeric_index_type my_size       = this->size();
  const numeric_index_type my_local_size = (last_local_idx - first_local_idx + 1);
  PetscErrorCode ierr=0;

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
    IS is;
    VecScatter scatter;

    // Create idx, idx[i] = i+first_local_idx;
    std::vector<PetscInt> idx(my_local_size);
    Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateLibMesh(this->comm().get(), my_local_size,
                           my_local_size ? &idx[0] : libmesh_nullptr, PETSC_USE_POINTER, &is);
    LIBMESH_CHKERR(ierr);

    ierr = VecScatterCreate(_vec,              is,
                            parallel_vec._vec, is,
                            &scatter);
    LIBMESH_CHKERR(ierr);

    ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,
                           INSERT_VALUES, SCATTER_FORWARD);
    LIBMESH_CHKERR(ierr);

    ierr = VecScatterEnd  (scatter, _vec, parallel_vec._vec,
                           INSERT_VALUES, SCATTER_FORWARD);
    LIBMESH_CHKERR(ierr);

    // Clean up
    ierr = LibMeshISDestroy (&is);
    LIBMESH_CHKERR(ierr);

    ierr = LibMeshVecScatterDestroy(&scatter);
    LIBMESH_CHKERR(ierr);
  }

  // localize like normal
  parallel_vec.close();
  parallel_vec.localize (*this, send_list);
  this->close();
}



template <typename T>
void PetscVector<T>::localize (std::vector<T> & v_local) const
{
  this->_restore_array();

  // This function must be run on all processors at once
  parallel_object_only();

  PetscErrorCode ierr=0;
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
                                         const processor_id_type pid) const
{
  this->_restore_array();

  PetscErrorCode ierr=0;
  const PetscInt n  = size();
  const PetscInt nl = local_size();
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
  else
    {
      if(pid == 0) // optimized version for localizing to 0
        {
          Vec vout;
          VecScatter ctx;

          ierr = VecScatterCreateToZero(_vec, &ctx, &vout);
          LIBMESH_CHKERR(ierr);

          ierr = VecScatterBegin(ctx, _vec, vout, INSERT_VALUES, SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);
          ierr = VecScatterEnd(ctx, _vec, vout, INSERT_VALUES, SCATTER_FORWARD);
          LIBMESH_CHKERR(ierr);

          if(processor_id() == 0)
            {
              v_local.resize(n);

              ierr = VecGetArray (vout, &values);
              LIBMESH_CHKERR(ierr);

              for (PetscInt i=0; i<n; i++)
                v_local[i] = static_cast<Real>(values[i]);

              ierr = VecRestoreArray (vout, &values);
              LIBMESH_CHKERR(ierr);
            }

          ierr = LibMeshVecScatterDestroy(&ctx);
          LIBMESH_CHKERR(ierr);
          ierr = LibMeshVecDestroy(&vout);
          LIBMESH_CHKERR(ierr);

        }
      else
        {
          v_local.resize(n);

          numeric_index_type ioff = this->first_local_index();
          std::vector<Real> local_values (n, 0.);

          {
            ierr = VecGetArray (_vec, &values);
            LIBMESH_CHKERR(ierr);

            for (PetscInt i=0; i<nl; i++)
              local_values[i+ioff] = static_cast<Real>(values[i]);

            ierr = VecRestoreArray (_vec, &values);
            LIBMESH_CHKERR(ierr);
          }


          MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
                      pid, this->comm().get());
        }
    }
}

#endif


// Full specialization for Complex datatypes
#ifdef LIBMESH_USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize_to_one (std::vector<Complex> & v_local,
                                            const processor_id_type pid) const
{
  this->_restore_array();

  PetscErrorCode ierr=0;
  const PetscInt n  = size();
  const PetscInt nl = local_size();
  PetscScalar * values;


  v_local.resize(n);


  for (PetscInt i=0; i<n; i++)
    v_local[i] = 0.;

  // only one processor
  if (n == nl)
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
      MPI_Reduce (&real_local_values[0], &real_v_local[0], n,
                  MPI_REAL, MPI_SUM,
                  pid, this->comm().get());

      MPI_Reduce (&imag_local_values[0], &imag_v_local[0], n,
                  MPI_REAL, MPI_SUM,
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
  this->_restore_array();

  PetscErrorCode ierr = 0;

  // Convert arguments to PetscVector*.
  const PetscVector<T> * vec1_petsc = cast_ptr<const PetscVector<T> *>(&vec1);
  const PetscVector<T> * vec2_petsc = cast_ptr<const PetscVector<T> *>(&vec2);

  // Call PETSc function.

  if(this->type() != GHOSTED)
    {
      ierr = VecPointwiseMult(this->vec(),
                              const_cast<PetscVector<T> *>(vec1_petsc)->vec(),
                              const_cast<PetscVector<T> *>(vec2_petsc)->vec());
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      Vec loc_vec;
      Vec v1_loc_vec;
      Vec v2_loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostGetLocalForm (const_cast<PetscVector<T> *>(vec1_petsc)->vec(),&v1_loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostGetLocalForm (const_cast<PetscVector<T> *>(vec2_petsc)->vec(),&v2_loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecPointwiseMult(loc_vec,v1_loc_vec,v2_loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecGhostRestoreLocalForm (const_cast<PetscVector<T> *>(vec1_petsc)->vec(),&v1_loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostRestoreLocalForm (const_cast<PetscVector<T> *>(vec2_petsc)->vec(),&v2_loc_vec);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
    }
}



template <typename T>
void PetscVector<T>::print_matlab (const std::string & name) const
{
  this->_restore_array();
  libmesh_assert (this->closed());

  PetscErrorCode ierr=0;
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (this->comm().get(),
                            &petsc_viewer);
  LIBMESH_CHKERR(ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.
   */
  if (name != "")
    {
      ierr = PetscViewerASCIIOpen( this->comm().get(),
                                   name.c_str(),
                                   &petsc_viewer);
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

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
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


  /**
   * Destroy the viewer.
   */
  ierr = LibMeshPetscViewerDestroy (&petsc_viewer);
  LIBMESH_CHKERR(ierr);
}





template <typename T>
void PetscVector<T>::create_subvector(NumericVector<T> & subvector,
                                      const std::vector<numeric_index_type> & rows) const
{
  this->_restore_array();

  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  PetscErrorCode ierr = 0;

  // Make sure the passed in subvector is really a PetscVector
  PetscVector<T> * petsc_subvector = cast_ptr<PetscVector<T> *>(&subvector);

  // If the petsc_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate PETSc routines to initialize it.
  if (!petsc_subvector->initialized())
    {
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
      LIBMESH_CHKERR(ierr);

      ierr = VecSetFromOptions (petsc_subvector->_vec);
      LIBMESH_CHKERR(ierr);

      // Mark the subvector as initialized
      petsc_subvector->_is_initialized = true;
    }
  else
    {
      petsc_subvector->_restore_array();
    }

  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<PetscInt> idx(rows.size());
  Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateLibMesh(this->comm().get(),
                         rows.size(),
                         numeric_petsc_cast(&rows[0]),
                         PETSC_USE_POINTER,
                         &parent_is);
  LIBMESH_CHKERR(ierr);

  ierr = ISCreateLibMesh(this->comm().get(),
                         rows.size(),
                         &idx[0],
                         PETSC_USE_POINTER,
                         &subvector_is);
  LIBMESH_CHKERR(ierr);

  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,
                          parent_is,
                          petsc_subvector->_vec,
                          subvector_is,
                          &scatter); LIBMESH_CHKERR(ierr);

  // Actually perform the scatter
  ierr = VecScatterBegin(scatter,
                         this->_vec,
                         petsc_subvector->_vec,
                         INSERT_VALUES,
                         SCATTER_FORWARD); LIBMESH_CHKERR(ierr);

  ierr = VecScatterEnd(scatter,
                       this->_vec,
                       petsc_subvector->_vec,
                       INSERT_VALUES,
                       SCATTER_FORWARD); LIBMESH_CHKERR(ierr);

  // Clean up
  ierr = LibMeshISDestroy(&parent_is);       LIBMESH_CHKERR(ierr);
  ierr = LibMeshISDestroy(&subvector_is);    LIBMESH_CHKERR(ierr);
  ierr = LibMeshVecScatterDestroy(&scatter); LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscVector<T>::_get_array(bool read_only) const
{
  libmesh_assert (this->initialized());

#ifdef LIBMESH_HAVE_CXX11_THREAD
  std::atomic_thread_fence(std::memory_order_acquire);
#else
  Threads::spin_mutex::scoped_lock lock(_petsc_vector_mutex);
#endif

  // If the values have already been retrieved and we're currently
  // trying to get a non-read only view (ie read/write) and the
  // values are currently read only... then we need to restore
  // the array first... and then retrieve it again.
  if (_array_is_present && !read_only && _values_read_only)
    _restore_array();

  // If we already have a read/write array - and we're trying
  // to get a read only array - let's just use the read write
  if (_array_is_present && read_only && !_values_read_only)
    _read_only_values = _values;

  if (!_array_is_present)
    {
#ifdef LIBMESH_HAVE_CXX11_THREAD
      std::lock_guard<std::mutex> lock(_petsc_vector_mutex);
#endif
      if (!_array_is_present)
        {
          PetscErrorCode ierr=0;
          if(this->type() != GHOSTED)
            {
#if PETSC_VERSION_LESS_THAN(3,2,0)
              // Vec{Get,Restore}ArrayRead were introduced in PETSc 3.2.0.  If you
              // have an older PETSc than that, we'll do an ugly
              // const_cast and call VecGetArray() instead.
              ierr = VecGetArray(_vec, const_cast<PetscScalar **>(&_values));
              _values_read_only = false;
#else
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
#endif
              LIBMESH_CHKERR(ierr);
              _local_size = this->local_size();
            }
          else
            {
              ierr = VecGhostGetLocalForm (_vec,&_local_form);
              LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,2,0)
              // Vec{Get,Restore}ArrayRead were introduced in PETSc 3.2.0.  If you
              // have an older PETSc than that, we'll do an ugly
              // const_cast and call VecGetArray() instead.
              ierr = VecGetArray(_local_form, const_cast<PetscScalar **>(&_values));
              _values_read_only = false;
#else
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
#endif
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
#ifdef LIBMESH_HAVE_CXX11_THREAD
          std::atomic_thread_fence(std::memory_order_release);
          _array_is_present.store(true, std::memory_order_relaxed);
#else
          _array_is_present = true;
#endif
        }
    }
}



template <typename T>
void PetscVector<T>::_restore_array() const
{
  if (_values_manually_retrieved)
    libmesh_error_msg("PetscVector values were manually retrieved but are being automatically restored!");

#ifdef LIBMESH_HAVE_CXX11_THREAD
  std::atomic_thread_fence(std::memory_order_acquire);
#else
  Threads::spin_mutex::scoped_lock lock(_petsc_vector_mutex);
#endif

  libmesh_assert (this->initialized());
  if (_array_is_present)
    {
#ifdef LIBMESH_HAVE_CXX11_THREAD
      std::lock_guard<std::mutex> lock(_petsc_vector_mutex);
#endif

      if (_array_is_present)
        {
          PetscErrorCode ierr=0;
          if(this->type() != GHOSTED)
            {
#if PETSC_VERSION_LESS_THAN(3,2,0)
              // Vec{Get,Restore}ArrayRead were introduced in PETSc 3.2.0.  If you
              // have an older PETSc than that, we'll do an ugly
              // const_cast and call VecRestoreArray() instead.
              ierr = VecRestoreArray (_vec, const_cast<PetscScalar **>(&_values));
#else
              if (_values_read_only)
                ierr = VecRestoreArrayRead (_vec, &_read_only_values);
              else
                ierr = VecRestoreArray (_vec, &_values);
#endif

              LIBMESH_CHKERR(ierr);
              _values = libmesh_nullptr;
            }
          else
            {
#if PETSC_VERSION_LESS_THAN(3,2,0)
              // Vec{Get,Restore}ArrayRead were introduced in PETSc 3.2.0.  If you
              // have an older PETSc than that, we'll do an ugly
              // const_cast and call VecRestoreArray() instead.
              ierr = VecRestoreArray (_local_form, const_cast<PetscScalar **>(&_values));
#else
              if (_values_read_only)
                ierr = VecRestoreArrayRead (_local_form, &_read_only_values);
              else
                ierr = VecRestoreArray (_local_form, &_values);
#endif
              LIBMESH_CHKERR(ierr);
              _values = libmesh_nullptr;
              ierr = VecGhostRestoreLocalForm (_vec,&_local_form);
              LIBMESH_CHKERR(ierr);
              _local_form = libmesh_nullptr;
              _local_size = 0;
            }
#ifdef LIBMESH_HAVE_CXX11_THREAD
          std::atomic_thread_fence(std::memory_order_release);
          _array_is_present.store(false, std::memory_order_relaxed);
#else
          _array_is_present = false;
#endif
        }
    }
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscVector<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC
