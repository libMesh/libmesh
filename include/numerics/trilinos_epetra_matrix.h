// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TRILINOS_EPETRA_MATRIX_H
#define LIBMESH_TRILINOS_EPETRA_MATRIX_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_TRILINOS

// Trilinos includes
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

// The EpetraExt interface is only needed for EpetraMatrix::add()
#ifdef LIBMESH_TRILINOS_HAVE_EPETRAEXT
#  include <EpetraExt_MatrixMatrix.h>
#endif

// Local includes
#include "libmesh/sparse_matrix.h"

// C++ includes
#include <algorithm>
#include <cstddef>

namespace libMesh
{

// Forward Declarations
template <typename T> class DenseMatrix;



/**
 * Epetra matrix. Provides a nice interface to the
 * Epetra data structures for parallel,
 * sparse matrices.
 *
 * \author Benjamin S. Kirk
 * \date 2008
 */
template <typename T>
class EpetraMatrix libmesh_final : public SparseMatrix<T>
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
  EpetraMatrix (const Parallel::Communicator & comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Constructor.  Creates a EpetraMatrix assuming you already
   * have a valid Epetra_FECrsMatrix object.  In this case, m is NOT destroyed
   * by the EpetraMatrix destructor when this object goes out of scope.
   * This allows ownership of m to remain with the original creator,
   * and to simply provide additional functionality with the EpetraMatrix.
   */
  EpetraMatrix (Epetra_FECrsMatrix * m,
                const Parallel::Communicator & comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  virtual ~EpetraMatrix ();

  /**
   * The \p EpetraMatrix needs the full sparsity pattern.
   */
  virtual bool need_full_sparsity_pattern () const libmesh_override
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will
   * tell the underlying matrix storage scheme how
   * to map the \f$ (i,j) \f$ elements.
   */
  virtual void update_sparsity_pattern (const SparsityPattern::Graph &) libmesh_override;

  /**
   * Initialize a Petsc matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * Optionally supports a block size, which indicates dense coupled blocks
   * for systems with multiple variables all of the same type.
   */
  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type nnz=30,
                     const numeric_index_type noz=10,
                     const numeric_index_type blocksize=1) libmesh_override;

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */
  virtual void init () libmesh_override;

  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor.
   */
  virtual void clear () libmesh_override;

  /**
   * Set all entries to 0. This method retains
   * sparsity structure.
   */
  virtual void zero () libmesh_override;

  /**
   * Call the Petsc assemble routines.
   * sends necessary messages to other
   * processors
   */
  virtual void close () const libmesh_override;

  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  virtual numeric_index_type m () const libmesh_override;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  virtual numeric_index_type n () const libmesh_override;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  virtual numeric_index_type row_start () const libmesh_override;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  virtual numeric_index_type row_stop () const libmesh_override;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  virtual void set (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  virtual void add (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  /**
   * Add the full matrix to the
   * Petsc matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & rows,
                           const std::vector<numeric_index_type> & cols) libmesh_override;

  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this:
   * \f$\texttt{this} = a*X + \texttt{this} \f$.
   * It is advisable to not only allocate appropriate memory with
   * \p init() , but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.  Note: \p X will
   * be closed, if not already done, before performing any work.
   */
  virtual void add (const T a, SparseMatrix<T> & X) libmesh_override;

  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation, and you
   * should always be careful where
   * you call this function.
   */
  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const libmesh_override;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  virtual Real l1_norm () const libmesh_override;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_infty=max_{all rows i}\sum_{all
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_infty \leq |M|_infty |v|_infty\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  virtual Real linfty_norm () const libmesh_override;

  /**
   * see if Petsc matrix has been closed
   * and fully assembled yet
   */
  virtual bool closed() const libmesh_override;

  /**
   * Print the contents of the matrix, by default to libMesh::out.
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const libmesh_override;

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T> & dest) const libmesh_override;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMatrix<T> & dest) const libmesh_override;

  /**
   * Swaps the raw PETSc matrix context pointers.
   */
  void swap (EpetraMatrix<T> &);

  /**
   * Returns the raw PETSc matrix context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling LibMeshMatDestroy()!
   */
  Epetra_FECrsMatrix * mat () { libmesh_assert(_mat); return _mat; }

  const Epetra_FECrsMatrix * mat () const { libmesh_assert(_mat); return _mat; }


private:

  /**
   * Actual Epetra datatype
   * to hold matrix entries
   */
  Epetra_FECrsMatrix * _mat;

  /**
   * Holds the distributed Map
   */
  Epetra_Map * _map;

  /**
   * Holds the sparsity pattern
   */
  Epetra_CrsGraph * _graph;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Mat object.
   */
  bool _destroy_mat_on_exit;

  /**
   * Epetra has no GetUseTranspose so we need to keep track of whether
   * we're transposed manually.
   */
  bool _use_transpose;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_TRILINOS
#endif // LIBMESH_TRILINOS_EPETRA_MATRIX_H
