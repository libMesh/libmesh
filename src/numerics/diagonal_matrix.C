// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/diagonal_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(const Parallel::Communicator & comm_in) : SparseMatrix<T>(comm_in)
{
  _diagonal = NumericVector<T>::build(comm_in);
}

template <typename T>
DiagonalMatrix<T> &
DiagonalMatrix<T>::operator=(const NumericVector<T> & vec)
{
  *_diagonal = vec;
  return *this;
}

template <typename T>
DiagonalMatrix<T> &
DiagonalMatrix<T>::operator=(NumericVector<T> && vec)
{
  // Don't get confused by the &&: vec is an lvalue reference; the && just
  // indicates that we are receiving an object that is safe to move from. Note
  // that we are not going to use std::move here because we do not have
  // (virtual) move assignment operations defined for NumericVector sub-classes
  _diagonal->swap(vec);
  return *this;
}

template <typename T>
void
DiagonalMatrix<T>::init(const numeric_index_type m,
                        const numeric_index_type /*n*/,
                        const numeric_index_type m_l,
                        const numeric_index_type /*n_l*/,
                        const numeric_index_type /*nnz*/,
                        const numeric_index_type /*noz*/,
                        const numeric_index_type /*blocksize*/)
{
  _diagonal->init(m, m_l);
}

template <typename T>
void
DiagonalMatrix<T>::init(const ParallelType type)
{
  libmesh_assert(this->_dof_map);

  _diagonal->init(this->_dof_map->n_dofs(),
                  this->_dof_map->n_dofs_on_processor(this->processor_id()),
                  /*fast=*/false,
                  type);
}

template <typename T>
void
DiagonalMatrix<T>::init(const NumericVector<T> & other, const bool fast)
{
  _diagonal->init(other, fast);
}

template <typename T>
void
DiagonalMatrix<T>::init(const DiagonalMatrix<T> & other, const bool fast)
{
  init(other.diagonal(), fast);
}

template <typename T>
void
DiagonalMatrix<T>::clear()
{
  _diagonal->clear();
}

template <typename T>
void
DiagonalMatrix<T>::zero()
{
  _diagonal->zero();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>> DiagonalMatrix<T>::zero_clone () const
{
  // Make empty copy with matching comm
  auto mat_copy = libmesh_make_unique<DiagonalMatrix<T>>(this->comm());

  // Initialize copy with our same nonzero structure, and explicitly
  // zero values using fast == false.
  mat_copy->init(*this, /*fast=*/false);

  // Work around an issue on older compilers.  We are able to simply
  // "return mat_copy;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(mat_copy.release());
}



template <typename T>
std::unique_ptr<SparseMatrix<T>> DiagonalMatrix<T>::clone () const
{
  // Make empty copy with matching comm
  auto mat_copy = libmesh_make_unique<DiagonalMatrix<T>>(this->comm());

  // Make copy of our diagonal
  auto diag_copy = _diagonal->clone();

  // Swap diag_copy with diagonal in mat_copy
  *mat_copy = std::move(*diag_copy);

  // Work around an issue on older compilers.  We are able to simply
  // "return mat_copy;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(mat_copy.release());
}

template <typename T>
void
DiagonalMatrix<T>::close()
{
  _diagonal->close();
}

template <typename T>
numeric_index_type
DiagonalMatrix<T>::m() const
{
  return _diagonal->size();
}

template <typename T>
numeric_index_type
DiagonalMatrix<T>::n() const
{
  return _diagonal->size();
}

template <typename T>
numeric_index_type
DiagonalMatrix<T>::row_start() const
{
  return _diagonal->first_local_index();
}

template <typename T>
numeric_index_type
DiagonalMatrix<T>::row_stop() const
{
  return _diagonal->last_local_index();
}

template <typename T>
void
DiagonalMatrix<T>::set(const numeric_index_type i, const numeric_index_type j, const T value)
{
  if (i == j)
    _diagonal->set(i, value);
}

template <typename T>
void
DiagonalMatrix<T>::add(const numeric_index_type i, const numeric_index_type j, const T value)
{
  if (i == j)
    _diagonal->add(i, value);
}

template <typename T>
void
DiagonalMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                              const std::vector<numeric_index_type> & rows,
                              const std::vector<numeric_index_type> & cols)
{
  auto m = dm.m();
  auto n = dm.n();

  for (decltype(m) i = 0; i < m; ++i)
    for (decltype(n) j = 0; j < n; ++j)
    {
      auto global_i = rows[i];
      auto global_j = cols[j];
      if (global_i == global_j)
        _diagonal->add(global_i, dm(i, j));
    }
}

template <typename T>
void
DiagonalMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                              const std::vector<numeric_index_type> & dof_indices)
{
  _diagonal->add_vector(dm.diagonal(), dof_indices);
}

template <typename T>
void
DiagonalMatrix<T>::add(const T a, const SparseMatrix<T> & X)
{
  auto x_diagonal = _diagonal->zero_clone();
  X.get_diagonal(*x_diagonal);
  _diagonal->add(a, *x_diagonal);
}

template <typename T>
T
DiagonalMatrix<T>::operator()(const numeric_index_type i, const numeric_index_type j) const
{
  if (i == j)
    return (*_diagonal)(i);
  else
    return 0;
}

template <typename T>
Real
DiagonalMatrix<T>::l1_norm() const
{
  return _diagonal->l1_norm();
}

template <typename T>
Real
DiagonalMatrix<T>::linfty_norm() const
{
  return _diagonal->linfty_norm();
}

template <typename T>
bool
DiagonalMatrix<T>::closed() const
{
  return _diagonal->closed();
}

template <typename T>
void
DiagonalMatrix<T>::print_personal(std::ostream & os) const
{
  _diagonal->print(os);
}

template <typename T>
void
DiagonalMatrix<T>::get_diagonal(NumericVector<T> & dest) const
{
  dest = *_diagonal;
}

template <typename T>
void
DiagonalMatrix<T>::get_transpose(SparseMatrix<T> & dest) const
{
  auto diagonal_dest = dynamic_cast<DiagonalMatrix<T> *>(&dest);
  if (diagonal_dest)
    *diagonal_dest = *_diagonal;
  else
    libmesh_error_msg("DenseMatrix<T>::get_transpose currently only accepts another DenseMatrix<T> "
                      "as its argument");
}

template <typename T>
void
DiagonalMatrix<T>::zero_rows(std::vector<numeric_index_type> & rows, T val/*=0*/)
{
  for (auto row : rows)
    _diagonal->set(row, val);
}

template <typename T>
const NumericVector<T> &
DiagonalMatrix<T>::diagonal() const
{
  return *_diagonal;
}

template class LIBMESH_EXPORT DiagonalMatrix<Number>;
}
