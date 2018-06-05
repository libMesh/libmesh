// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LASPACK_MATRIX_H
#define LIBMESH_LASPACK_MATRIX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_LASPACK

// Local includes
#include "libmesh/sparse_matrix.h"

// Laspack includes
#include <qmatrix.h>

// C++ includes
#include <algorithm>
#include <cstddef>

namespace libMesh
{

// Forward declarations
template <typename T> class DenseMatrix;
template <typename T> class LaspackVector;
template <typename T> class LaspackLinearSolver;

/**
 * The LaspackMatrix class wraps a QMatrix object from the Laspack
 * library. Currently Laspack only supports real datatypes, so this
 * class is a full specialization of \p SparseMatrix<T> with \p T = \p
 * Real. All overridden virtual functions are documented in
 * sparse_matrix.h.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template <typename T>
class LaspackMatrix libmesh_final : public SparseMatrix<T>
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
   * You have to initialize the matrix before usage with \p init(...).
   */
  LaspackMatrix (const Parallel::Communicator & comm);

  /**
   * Destructor. Free all memory, but do not release the memory of the
   * sparsity structure.
   */
  ~LaspackMatrix ();

  /**
   * The \p LaspackMatrix needs the full sparsity pattern.
   */
  virtual bool need_full_sparsity_pattern() const override
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will tell the
   * underlying matrix storage scheme how to map the \f$ (i,j) \f$
   * elements.
   */
  virtual void update_sparsity_pattern (const SparsityPattern::Graph &) override;

  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type nnz=30,
                     const numeric_index_type noz=10,
                     const numeric_index_type blocksize=1) override;

  virtual void init () override;

  virtual void clear () override;

  virtual void zero () override;

  virtual void close () override { this->_closed = true; }

  virtual numeric_index_type m () const override;

  virtual numeric_index_type n () const override;

  virtual numeric_index_type row_start () const override;

  virtual numeric_index_type row_stop () const override;

  virtual void set (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) override;

  virtual void add (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & rows,
                           const std::vector<numeric_index_type> & cols) override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) override;

  /**
   * Compute A += a*X for scalar \p a, matrix \p X.
   *
   * \note LASPACK does not provide a true \p axpy for matrices,
   * so a hand-coded version with hopefully acceptable performance
   * is provided.
   */
  virtual void add (const T a, SparseMatrix<T> & X) override;

  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const override;

  virtual Real l1_norm () const override { libmesh_not_implemented(); return 0.; }

  virtual Real linfty_norm () const override { libmesh_not_implemented(); return 0.; }

  virtual bool closed() const override { return _closed; }

  virtual void print_personal(std::ostream & os=libMesh::out) const override { this->print(os); }

  virtual void get_diagonal (NumericVector<T> & dest) const override;

  virtual void get_transpose (SparseMatrix<T> & dest) const override;

private:

  /**
   * \returns The position in the compressed row
   * storage scheme of the \f$ (i,j) \f$ element.
   */
  numeric_index_type pos (const numeric_index_type i,
                          const numeric_index_type j) const;

  /**
   *  The Laspack sparse matrix pointer.
   */
  QMatrix _QMat;

  /**
   * The compressed row indices.
   */
  std::vector<numeric_index_type> _csr;

  /**
   * The start of each row in the compressed
   * row index data structure.
   */
  std::vector<std::vector<numeric_index_type>::const_iterator> _row_start;

  /**
   * Flag indicating if the matrix has been closed yet.
   */
  bool _closed;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackVector<T>;
  friend class LaspackLinearSolver<T>;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifdef LIBMESH_LASPACK_MATRIX_H
