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



#ifndef LIBMESH_TRILINOS_EPETRA_MATRIX_H
#define LIBMESH_TRILINOS_EPETRA_MATRIX_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

// Trilinos includes
#include "libmesh/ignore_warnings.h"
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

// The EpetraExt interface is only needed for EpetraMatrix::add()
#ifdef LIBMESH_TRILINOS_HAVE_EPETRAEXT
#  include <EpetraExt_MatrixMatrix.h>
#endif
#include "libmesh/restore_warnings.h"

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
 * This class provides a nice interface to the Epetra data structures
 * for parallel, sparse matrices. All overridden virtual functions are
 * documented in sparse_matrix.h.
 *
 * \author Benjamin S. Kirk
 * \date 2008
 */
template <typename T>
class EpetraMatrix libmesh_final : public SparseMatrix<T>
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
  EpetraMatrix (const Parallel::Communicator & comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Constructor.  Creates a EpetraMatrix assuming you already have a
   * valid Epetra_FECrsMatrix object.  In this case, m is NOT
   * destroyed by the EpetraMatrix destructor when this object goes
   * out of scope.  This allows ownership of m to remain with the
   * original creator, and to simply provide additional functionality
   * with the EpetraMatrix.
   */
  EpetraMatrix (Epetra_FECrsMatrix * m,
                const Parallel::Communicator & comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor. Free all memory, but do not release the memory of the
   * sparsity structure.
   */
  virtual ~EpetraMatrix ();

  /**
   * The \p EpetraMatrix needs the full sparsity pattern.
   */
  virtual bool need_full_sparsity_pattern () const libmesh_override
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will tell the
   * underlying matrix storage scheme how to map the \f$ (i,j) \f$
   * elements.
   */
  virtual void update_sparsity_pattern (const SparsityPattern::Graph &) libmesh_override;

  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type nnz=30,
                     const numeric_index_type noz=10,
                     const numeric_index_type blocksize=1) libmesh_override;

  virtual void init () libmesh_override;

  virtual void clear () libmesh_override;

  virtual void zero () libmesh_override;

  virtual void close () libmesh_override;

  virtual numeric_index_type m () const libmesh_override;

  virtual numeric_index_type n () const libmesh_override;

  virtual numeric_index_type row_start () const libmesh_override;

  virtual numeric_index_type row_stop () const libmesh_override;

  virtual void set (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  virtual void add (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & rows,
                           const std::vector<numeric_index_type> & cols) libmesh_override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  /**
   * Compute A += a*X for scalar \p a, matrix \p X.
   *
   * \note It is advisable to not only allocate appropriate memory
   * with \p init(), but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.
   *
   * \note \p X will be closed, if not already done, before performing
   * any work.
   */
  virtual void add (const T a, SparseMatrix<T> & X) libmesh_override;

  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const libmesh_override;

  virtual Real l1_norm () const libmesh_override;

  virtual Real linfty_norm () const libmesh_override;

  virtual bool closed() const libmesh_override;

  virtual void print_personal(std::ostream & os=libMesh::out) const libmesh_override;

  virtual void get_diagonal (NumericVector<T> & dest) const libmesh_override;

  virtual void get_transpose (SparseMatrix<T> & dest) const libmesh_override;

  /**
   * Swaps the internal data pointers, no actual values are swapped.
   */
  void swap (EpetraMatrix<T> &);

  /**
   * \returns The raw Epetra_FECrsMatrix pointer.
   *
   * \note This is generally not required in user-level code.
   *
   * \note Don't do anything crazy like deleting the pointer, or very
   * bad things will likely happen!
   */
  Epetra_FECrsMatrix * mat () { libmesh_assert(_mat); return _mat; }

  const Epetra_FECrsMatrix * mat () const { libmesh_assert(_mat); return _mat; }


private:

  /**
   * Actual Epetra datatype to hold matrix entries.
   */
  Epetra_FECrsMatrix * _mat;

  /**
   * Holds the distributed Map.
   */
  Epetra_Map * _map;

  /**
   * Holds the sparsity pattern.
   */
  Epetra_CrsGraph * _graph;

  /**
   * This boolean value should only be set to false for the
   * constructor which takes an Epetra_FECrsMatrix object.
   */
  bool _destroy_mat_on_exit;

  /**
   * Epetra has no GetUseTranspose so we need to keep track of whether
   * we're transposed manually.
   */
  bool _use_transpose;
};

} // namespace libMesh

#endif // LIBMESH_TRILINOS_HAVE_EPETRA
#endif // LIBMESH_TRILINOS_EPETRA_MATRIX_H
