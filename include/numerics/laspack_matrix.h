// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
 * Real.
 *
 * @author Benjamin S. Kirk, 2003
 */

template <typename T>
class LaspackMatrix : public SparseMatrix<T>
{

public:
  /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */
  LaspackMatrix (const Parallel::Communicator &comm = libMesh::CommWorld);

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~LaspackMatrix ();

  /**
   * The \p LaspackMatrix needs the full sparsity pattern.
   */
  bool need_full_sparsity_pattern() const
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will
   * tell the underlying matrix storage scheme how
   * to map the \f$ (i,j) \f$ elements.
   */
  void update_sparsity_pattern (const SparsityPattern::Graph &);

  /**
   * Initialize a Laspack matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   */
  void init (const numeric_index_type m,
	     const numeric_index_type n,
	     const numeric_index_type m_l,
	     const numeric_index_type n_l,
	     const numeric_index_type nnz=30,
	     const numeric_index_type noz=10);

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */
  void init ();

  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor.
   */
  void clear ();

  /**
   * Set all entries to 0.
   */
  void zero ();

  /**
   * Close the matrix.  Dummy routine.  After calling
   * this method \p closed() is true and the matrix can
   * be used in computations.
   */
  void close () const { const_cast<LaspackMatrix<T>*>(this)->_closed = true; }

  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  numeric_index_type m () const;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  numeric_index_type n () const;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  numeric_index_type row_start () const;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  numeric_index_type row_stop () const;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  void set (const numeric_index_type i,
	    const numeric_index_type j,
	    const T value);

  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  void add (const numeric_index_type i,
	    const numeric_index_type j,
	    const T value);

  /**
   * Add the full matrix to the
   * Laspack matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */

  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<numeric_index_type> &rows,
		   const std::vector<numeric_index_type> &cols);

  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<numeric_index_type> &dof_indices);

  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this: \f$\texttt{this} += a*X \f$.
   * \p LASPACK does not provide a true \p axpy for matrices,
   * so a hand-coded version with hopefully acceptable performance
   * is provided.
   */
  void add (const T a, SparseMatrix<T> &X);

  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation, and you
   * should always be careful where
   * you call this function.
   */
  T operator () (const numeric_index_type i,
		 const numeric_index_type j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   */
  Real l1_norm () const { libmesh_error(); return 0.; }

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_\infty=max_{all rows i}\sum_{all
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_\infty \leq |M|_\infty |v|_\infty\f$.
   */
  Real linfty_norm () const { libmesh_error(); return 0.; }

  /**
   * see if Laspack matrix has been closed
   * and fully assembled yet
   */
  bool closed() const { return _closed; }

  /**
   * Print the contents of the matrix, by default to libMesh::out.
   * Currently identical to \p print().
   */
  void print_personal(std::ostream& os=libMesh::out) const { this->print(os); }

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T>& dest) const;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMatrix<T>& dest) const;

private:

  /**
   * @returns the position in the compressed row
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
