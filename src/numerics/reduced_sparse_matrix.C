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

#include "libmesh/reduced_sparse_matrix.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/reduced_dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dof_object.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{

template <typename T>
ReducedSparseMatrix<T>::ReducedSparseMatrix(const Parallel::Communicator & comm_in,
                                            const ReducedDofMap & reduced_dof_map)
  : PetscMatrixBase<T>(comm_in), _reduced_dof_map(reduced_dof_map)
{
  // The wrapped matrix owns the Mat
  this->_destroy_mat_on_exit = false;
}

template <typename T>
ReducedSparseMatrix<T>::~ReducedSparseMatrix() = default;

template <typename T>
void ReducedSparseMatrix<T>::adopt_reduced_matrix()
{
  this->_mat = _reduced->mat();
  this->_is_initialized = true;
  // Solvers recover this object from the Mat, so it, not the wrapped matrix, is the context
  this->set_context();
}

template <typename T>
void ReducedSparseMatrix<T>::init(ParallelType)
{
  libmesh_assert(_reduced_dof_map.initialized());

  const numeric_index_type n = _reduced_dof_map.n_dofs();
  const numeric_index_type n_local = _reduced_dof_map.n_local_dofs();

  _reduced = std::make_unique<PetscMatrix<T>>(this->comm());
  _reduced->init(
      n, n, n_local, n_local, _reduced_dof_map.reduced_n_nz(), _reduced_dof_map.reduced_n_oz());
  this->adopt_reduced_matrix();
}

template <typename T>
void ReducedSparseMatrix<T>::init(const numeric_index_type m,
                                  const numeric_index_type n,
                                  const numeric_index_type m_l,
                                  const numeric_index_type n_l,
                                  const numeric_index_type nnz,
                                  const numeric_index_type noz,
                                  const numeric_index_type blocksize)
{
  _reduced = std::make_unique<PetscMatrix<T>>(this->comm());
  _reduced->init(m, n, m_l, n_l, nnz, noz, blocksize);
  this->adopt_reduced_matrix();
}

template <typename T>
void ReducedSparseMatrix<T>::clear() noexcept
{
  _reduced.reset();
  this->_mat = nullptr;
  this->_is_initialized = false;
}

template <typename T>
PetscMatrix<T> & ReducedSparseMatrix<T>::reduced_matrix()
{
  libmesh_assert(_reduced);
  return *_reduced;
}

template <typename T>
SparseMatrix<T> & ReducedSparseMatrix<T>::operator=(const SparseMatrix<T> & v)
{
  const auto & other = cast_ref<const ReducedSparseMatrix<T> &>(v);
  libmesh_assert(&other._reduced_dof_map == &_reduced_dof_map);
  *_reduced = *other._reduced;
  return *this;
}

template <typename T>
void ReducedSparseMatrix<T>::to_reduced(const std::vector<numeric_index_type> & full,
                                        std::vector<PetscInt> & reduced) const
{
  reduced.resize(full.size());
  for (const auto i : index_range(full))
    {
      const dof_id_type r = _reduced_dof_map.find_reduced_dof(full[i]);
      reduced[i] = (r == DofObject::invalid_id) ? PetscInt(-1) : static_cast<PetscInt>(r);
    }
}

template <typename T>
void ReducedSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                        const std::vector<numeric_index_type> & rows,
                                        const std::vector<numeric_index_type> & cols)
{
  libmesh_assert(this->initialized());
  libmesh_assert_equal_to(rows.size(), dm.m());
  libmesh_assert_equal_to(cols.size(), dm.n());

  if (rows.empty() || cols.empty())
    return;

  this->to_reduced(rows, _rows);
  this->to_reduced(cols, _cols);
  LibmeshPetscCall(MatSetValues(this->_mat,
                                cast_int<PetscInt>(_rows.size()),
                                _rows.data(),
                                cast_int<PetscInt>(_cols.size()),
                                _cols.data(),
                                pPS(const_cast<T *>(dm.get_values().data())),
                                ADD_VALUES));
}

template <typename T>
void ReducedSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                        const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix(dm, dof_indices, dof_indices);
}

template <typename T>
void ReducedSparseMatrix<T>::set(const numeric_index_type i,
                                 const numeric_index_type j,
                                 const T value)
{
  const dof_id_type ri = _reduced_dof_map.find_reduced_dof(i);
  const dof_id_type rj = _reduced_dof_map.find_reduced_dof(j);
  if (ri != DofObject::invalid_id && rj != DofObject::invalid_id)
    _reduced->set(ri, rj, value);
}

template <typename T>
void ReducedSparseMatrix<T>::add(const numeric_index_type i,
                                 const numeric_index_type j,
                                 const T value)
{
  const dof_id_type ri = _reduced_dof_map.find_reduced_dof(i);
  const dof_id_type rj = _reduced_dof_map.find_reduced_dof(j);
  if (ri != DofObject::invalid_id && rj != DofObject::invalid_id)
    _reduced->add(ri, rj, value);
}

template <typename T>
void ReducedSparseMatrix<T>::zero()
{
  _reduced->zero();
}

template <typename T>
void ReducedSparseMatrix<T>::close()
{
  _reduced->close();
}

template <typename T>
void ReducedSparseMatrix<T>::flush()
{
  _reduced->flush();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>> ReducedSparseMatrix<T>::zero_clone() const
{
  auto ret = std::make_unique<ReducedSparseMatrix<T>>(this->comm(), _reduced_dof_map);
  ret->_reduced.reset(cast_ptr<PetscMatrix<T> *>(_reduced->zero_clone().release()));
  ret->adopt_reduced_matrix();
  return ret;
}

template <typename T>
std::unique_ptr<SparseMatrix<T>> ReducedSparseMatrix<T>::clone() const
{
  auto ret = std::make_unique<ReducedSparseMatrix<T>>(this->comm(), _reduced_dof_map);
  ret->_reduced.reset(cast_ptr<PetscMatrix<T> *>(_reduced->clone().release()));
  ret->adopt_reduced_matrix();
  return ret;
}

template <typename T>
void ReducedSparseMatrix<T>::add(const T a, const SparseMatrix<T> & X)
{
  const auto & other = cast_ref<const ReducedSparseMatrix<T> &>(X);
  _reduced->add(a, *other._reduced);
}

template <typename T>
T ReducedSparseMatrix<T>::operator()(const numeric_index_type i, const numeric_index_type j) const
{
  return (*_reduced)(i, j);
}

template <typename T>
Real ReducedSparseMatrix<T>::l1_norm() const
{
  return _reduced->l1_norm();
}

template <typename T>
Real ReducedSparseMatrix<T>::linfty_norm() const
{
  return _reduced->linfty_norm();
}

template <typename T>
void ReducedSparseMatrix<T>::print_personal(std::ostream & os) const
{
  _reduced->print_personal(os);
}

template <typename T>
void ReducedSparseMatrix<T>::get_diagonal(NumericVector<T> & dest) const
{
  _reduced->get_diagonal(dest);
}

template <typename T>
void ReducedSparseMatrix<T>::get_transpose(SparseMatrix<T> & dest) const
{
  auto & other = cast_ref<ReducedSparseMatrix<T> &>(dest);
  libmesh_assert(&other._reduced_dof_map == &_reduced_dof_map);
  if (!other._reduced)
    other._reduced = std::make_unique<PetscMatrix<T>>(this->comm());
  // The transpose replaces the destination's Mat, which the destination must then adopt
  _reduced->get_transpose(*other._reduced);
  other.adopt_reduced_matrix();
}

template <typename T>
void ReducedSparseMatrix<T>::get_row(numeric_index_type i,
                                     std::vector<numeric_index_type> & indices,
                                     std::vector<T> & values) const
{
  _reduced->get_row(i, indices, values);
}

template <typename T>
void ReducedSparseMatrix<T>::scale(const T scale)
{
  _reduced->scale(scale);
}

template class LIBMESH_EXPORT ReducedSparseMatrix<Number>;

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
