// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_DIAGONAL_MATRIX_H
#define LIBMESH_DIAGONAL_MATRIX_H

#include "libmesh/id_types.h"
#include "libmesh/sparse_matrix.h"

#include <memory>

namespace libMesh
{
template <typename>
class NumericVector;
namespace Parallel
{
class Communicator;
}

/**
 * Diagonal matrix class whose underlying storage is a vector
 *
 * \author Alexander D. Lindsay
 * \date 2019
 */
template <typename T>
class DiagonalMatrix : public SparseMatrix<T>
{
public:
  /**
   * Constructor; initializes the matrix to be empty, without any
   * structure, i.e.  the matrix is not usable at all. This
   * constructor is therefore only useful for matrices which are
   * members of a class. All other matrices should be created at a
   * point in the data flow where all necessary information is
   * available.
   *
   * You have to initialize the matrix before usage with
   * \p init(...).
   */
  explicit DiagonalMatrix(const Parallel::Communicator & comm);

  /**
   * unique pointers can be moved but not copied
   */
  DiagonalMatrix(DiagonalMatrix &&) = default;
  DiagonalMatrix & operator=(DiagonalMatrix &&) = default;

  /**
   * Copy contents from vec into underlying diagonal storage
   */
  DiagonalMatrix & operator=(const NumericVector<T> & vec);

  /**
   * Move contents from vec into underlying diagonal storage
   */
  DiagonalMatrix & operator=(NumericVector<T> && vec);


  void init(const numeric_index_type m,
            const numeric_index_type n,
            const numeric_index_type m_l,
            const numeric_index_type n_l,
            const numeric_index_type nnz = 30,
            const numeric_index_type noz = 10,
            const numeric_index_type blocksize = 1) override;

  void init() override;

  void clear() override;

  void zero() override;

  void close() override;

  numeric_index_type m() const override;

  numeric_index_type n() const override;

  numeric_index_type row_start() const override;

  numeric_index_type row_stop() const override;

  void set(const numeric_index_type i, const numeric_index_type j, const T value) override;

  void add(const numeric_index_type i, const numeric_index_type j, const T value) override;

  void add_matrix(const DenseMatrix<T> & dm,
                  const std::vector<numeric_index_type> & rows,
                  const std::vector<numeric_index_type> & cols) override;

  void add_matrix(const DenseMatrix<T> & dm,
                  const std::vector<numeric_index_type> & dof_indices) override;

  void add(const T a, const SparseMatrix<T> & X) override;

  T operator()(const numeric_index_type i, const numeric_index_type j) const override;

  Real l1_norm() const override;

  Real linfty_norm() const override;

  bool closed() const override;

  void print_personal(std::ostream & os = libMesh::out) const override;

  void get_diagonal(NumericVector<T> & dest) const override;

  void get_transpose(SparseMatrix<T> & dest) const override;

  void zero_rows(std::vector<numeric_index_type> & rows, T val = 0) override;

protected:
  /// Underlying diagonal matrix storage
  std::unique_ptr<NumericVector<T>> _diagonal;
};
} // namespace libMesh

#endif // LIBMESH_DIAGONAL_MATRIX_H
