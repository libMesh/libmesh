// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/lumped_mass_matrix.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/numeric_vector.h"
#include <cmath>

namespace libMesh
{

template <typename T>
LumpedMassMatrix<T>::LumpedMassMatrix(const Parallel::Communicator & comm_in)
  : DiagonalMatrix<T>(comm_in)
{
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
LumpedMassMatrix<T>::zero_clone() const
{
  // Make empty copy with matching comm
  auto mat_copy = libmesh_make_unique<LumpedMassMatrix<T>>(this->comm());

  // Initialize copy with our same nonzero structure, and explicitly
  // zero values using fast == false.
  mat_copy->init(*this, /*fast=*/false);

  // Work around an issue on older compilers.  We are able to simply
  // "return mat_copy;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(mat_copy.release());
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
LumpedMassMatrix<T>::clone() const
{
  // Make empty copy with matching comm
  auto mat_copy = libmesh_make_unique<LumpedMassMatrix<T>>(this->comm());

  // Make copy of our diagonal
  auto diag_copy = this->_diagonal->clone();

  // Swap diag_copy with diagonal in mat_copy
  *mat_copy = std::move(*diag_copy);

  // Work around an issue on older compilers.  We are able to simply
  // "return mat_copy;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(mat_copy.release());
}

template <typename T>
void
LumpedMassMatrix<T>::set(const numeric_index_type i,
                         const numeric_index_type libmesh_dbg_var(j),
                         const T value)
{
  libmesh_assert_msg(i == j, "Set in a lumped mass matrix really only makes sense for i == j");
  this->_diagonal->set(i, std::abs(value));
}

template <typename T>
void
LumpedMassMatrix<T>::add(const numeric_index_type i, const numeric_index_type /*j*/, const T value)
{
  this->_diagonal->add(i, std::abs(value));
}

template <typename T>
void
LumpedMassMatrix<T>::add(const T a, const SparseMatrix<T> & X)
{
  if (dynamic_cast<const LumpedMassMatrix<T> *>(&X))
  {
    auto x_diagonal = this->_diagonal->zero_clone();
    X.get_diagonal(*x_diagonal);
    this->_diagonal->add(a, *x_diagonal);
  }
  else
    libmesh_error_msg("Unsupported matrix type passed to LumpedMassMatrix::add");
}

template <typename T>
LumpedMassMatrix<T> &
LumpedMassMatrix<T>::operator=(const NumericVector<T> & vec)
{
  *this->_diagonal = vec;
  return *this;
}

template <typename T>
LumpedMassMatrix<T> &
LumpedMassMatrix<T>::operator=(NumericVector<T> && vec)
{
  this->_diagonal->swap(vec);
  return *this;
}

template class LumpedMassMatrix<Number>;

} // namespace libMesh
