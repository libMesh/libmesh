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



#ifndef LIBMESH_PETSC_MATRIX_H
#define LIBMESH_PETSC_MATRIX_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/parallel.h"

// C++ includes
#include <algorithm>

// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMatrix methods
#undef semiparallel_only
#ifndef NDEBUG
#include <cstring>

#define semiparallel_only() do { if (this->initialized()) { const char * mytype; \
      MatGetType(_mat,&mytype);                                         \
      if (!strcmp(mytype, MATSEQAIJ))                                   \
        parallel_object_only(); } } while (0)
#else
#define semiparallel_only()
#endif


// Petsc include files.
#include <petscmat.h>



namespace libMesh
{

// Forward Declarations
template <typename T> class DenseMatrix;


/**
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
template <typename T>
class PetscMatrix libmesh_final : public SparseMatrix<T>
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
  explicit
  PetscMatrix (const Parallel::Communicator & comm_in
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Constructor.  Creates a PetscMatrix assuming you already
   * have a valid Mat object.  In this case, m is NOT destroyed
   * by the PetscMatrix destructor when this object goes out of scope.
   * This allows ownership of m to remain with the original creator,
   * and to simply provide additional functionality with the PetscMatrix.
   */
  explicit
  PetscMatrix (Mat m,
               const Parallel::Communicator & comm_in
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~PetscMatrix ();

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
   * Initialize a Petsc matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row.
   * \p noz is the number of off-processor
   * nonzeros per row.
   * Optionally supports a block size, which indicates dense coupled blocks
   * for systems with multiple variables all of the same type.
   */
  void init (const numeric_index_type m,
             const numeric_index_type n,
             const numeric_index_type m_l,
             const numeric_index_type n_l,
             const std::vector<numeric_index_type> & n_nz,
             const std::vector<numeric_index_type> & n_oz,
             const numeric_index_type blocksize=1);

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */
  virtual void init () libmesh_override;

  /**
   * Update the sparsity pattern based on \p dof_map, and set the matrix
   * to zero. This is useful in cases where the sparsity pattern changes
   * during a computation.
   */
  void update_preallocation_and_zero();

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
   * Set all row entries to 0 then puts diag_value in the diagonal entry
   */
  virtual void zero_rows (std::vector<numeric_index_type> & rows, T diag_value = 0.0) libmesh_override;

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
   * Same as \p add_matrix, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  /**
   * Add the full matrix \p dm to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time.  The matrix is assumed blocked, and \p brow, \p bcol
   * correspond to the *block* row, columm indices.
   */
  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & brows,
                                 const std::vector<numeric_index_type> & bcols) libmesh_override;

  /**
   * Same as \p add_block_matrix , but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & dof_indices) libmesh_override
  { this->add_block_matrix (dm, dof_indices, dof_indices); }

  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this:
   * \f$\texttt{this} = a*X + \texttt{this} \f$.
   * Use this with caution, the sparse matrices need to have the
   * same nonzero pattern, otherwise \p PETSc will crash!
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
   * Print the contents of the matrix to the screen
   * with the PETSc viewer.  This function only allows
   * printing to standard out, this is because we have
   * limited ourselves to one PETSc implementation for
   * writing.
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const libmesh_override;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string & name = "") const libmesh_override;

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
  void swap (PetscMatrix<T> &);

  /**
   * Returns the raw PETSc matrix context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling LibMeshMatDestroy()!
   */
  Mat mat () { libmesh_assert (_mat); return _mat; }

protected:

  /**
   * This function either creates or re-initializes
   * a matrix called "submatrix" which is defined
   * by the row and column indices given in the "rows" and "cols" entries.
   * This function is implemented in terms of the MatGetSubMatrix()
   * routine of PETSc.  The boolean reuse_submatrix parameter determines
   * whether or not PETSc will treat "submatrix" as one which has already
   * been used (had memory allocated) or as a new matrix.
   */
  virtual void _get_submatrix(SparseMatrix<T> & submatrix,
                              const std::vector<numeric_index_type> & rows,
                              const std::vector<numeric_index_type> & cols,
                              const bool reuse_submatrix) const libmesh_override;

private:

  /**
   * Petsc matrix datatype to store values
   */
  Mat _mat;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Mat object.
   */
  bool _destroy_mat_on_exit;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_MATRIX_H
