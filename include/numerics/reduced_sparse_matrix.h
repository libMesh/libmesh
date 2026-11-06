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

#ifndef LIBMESH_REDUCED_SPARSE_MATRIX_H
#define LIBMESH_REDUCED_SPARSE_MATRIX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_matrix_base.h"

#include <memory>
#include <vector>

namespace libMesh
{
class ReducedDofMap;
template <typename T>
class PetscMatrix;

/**
 * A matrix on the reduced space of a \p ReducedDofMap, assembled with full degree of freedom
 * indices.
 *
 * Insertion (\p set(), \p add(), \p add_matrix()) takes indices of the full \p DofMap, because
 * element assembly produces them. Each index is mapped to its reduced index, and entries in rows
 * or columns that the \p ReducedDofMap does not retain are dropped. Dropping is exact when the
 * element contributions already carry the reduction, e.g. \f$ C^T K C \f$ for eliminated
 * constraints. Every other operation uses reduced indices: sizes, row and column ranges, entry
 * access, and the PETSc \p Mat returned by \p mat(), which is what solvers see.
 */
template <typename T>
class ReducedSparseMatrix : public PetscMatrixBase<T>
{
public:
  ReducedSparseMatrix(const Parallel::Communicator & comm_in, const ReducedDofMap & reduced_dof_map);

  virtual ~ReducedSparseMatrix();

  ReducedSparseMatrix(ReducedSparseMatrix &&) = delete;
  ReducedSparseMatrix(const ReducedSparseMatrix &) = delete;
  ReducedSparseMatrix & operator=(ReducedSparseMatrix &&) = delete;
  ReducedSparseMatrix & operator=(const ReducedSparseMatrix &) = delete;

  virtual SparseMatrix<T> & operator=(const SparseMatrix<T> & v) override;

  /**
   * Allocates the reduced matrix from the \p ReducedDofMap sparsity pattern
   */
  virtual void init(ParallelType type = PARALLEL) override;

  /**
   * Arguments are reduced sizes
   */
  virtual void init(const numeric_index_type m,
                    const numeric_index_type n,
                    const numeric_index_type m_l,
                    const numeric_index_type n_l,
                    const numeric_index_type nnz = 30,
                    const numeric_index_type noz = 10,
                    const numeric_index_type blocksize = 1) override;

  virtual bool require_sparsity_pattern() const override { return false; }

  virtual void clear() noexcept override;

  virtual void zero() override;

  virtual std::unique_ptr<SparseMatrix<T>> zero_clone() const override;

  virtual std::unique_ptr<SparseMatrix<T>> clone() const override;

  virtual void close() override;

  virtual void flush() override;

  /**
   * \p i and \p j are full indices
   */
  virtual void set(const numeric_index_type i, const numeric_index_type j, const T value) override;

  /**
   * \p i and \p j are full indices
   */
  virtual void add(const numeric_index_type i, const numeric_index_type j, const T value) override;

  /**
   * \p rows and \p cols are full indices
   */
  virtual void add_matrix(const DenseMatrix<T> & dm,
                          const std::vector<numeric_index_type> & rows,
                          const std::vector<numeric_index_type> & cols) override;

  /**
   * \p dof_indices are full indices
   */
  virtual void add_matrix(const DenseMatrix<T> & dm,
                          const std::vector<numeric_index_type> & dof_indices) override;

  virtual void add(const T a, const SparseMatrix<T> & X) override;

  /**
   * \p i and \p j are reduced indices
   */
  virtual T operator()(const numeric_index_type i, const numeric_index_type j) const override;

  virtual Real l1_norm() const override;

  virtual Real linfty_norm() const override;

  virtual void print_personal(std::ostream & os = libMesh::out) const override;

  /**
   * Fills the reduced-size vector \p dest with the diagonal
   */
  virtual void get_diagonal(NumericVector<T> & dest) const override;

  virtual void get_transpose(SparseMatrix<T> & dest) const override;

  /**
   * \p i and the returned \p indices are reduced indices
   */
  virtual void get_row(numeric_index_type i,
                       std::vector<numeric_index_type> & indices,
                       std::vector<T> & values) const override;

  virtual void scale(const T scale) override;

  /**
   * \returns The reduced degree of freedom map
   */
  const ReducedDofMap & reduced_dof_map() const { return _reduced_dof_map; }

  /**
   * \returns The wrapped matrix on the reduced space
   */
  PetscMatrix<T> & reduced_matrix();

private:
  /**
   * Maps the full indices \p full to reduced indices in \p reduced, with \p -1 for unretained
   * ones, which PETSc ignores on insertion
   */
  void to_reduced(const std::vector<numeric_index_type> & full,
                  std::vector<PetscInt> & reduced) const;

  /**
   * Wraps the reduced matrix's \p Mat and sets it as this matrix's context
   */
  void adopt_reduced_matrix();

  const ReducedDofMap & _reduced_dof_map;

  /// The matrix on the reduced space
  std::unique_ptr<PetscMatrix<T>> _reduced;

  /// Work vectors of reduced row and column indices
  mutable std::vector<PetscInt> _rows, _cols;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_REDUCED_SPARSE_MATRIX_H
