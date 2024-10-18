// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_solver_exception.h"

// C++ includes
#include <algorithm>
#ifdef LIBMESH_HAVE_CXX11_THREAD
#include <atomic>
#include <mutex>
#endif



namespace libMesh
{

// Forward Declarations
template <typename T> class DenseMatrix;

enum PetscMatrixType : int {
                 AIJ=0,
                 HYPRE};

/**
 * This class provides a nice interface to the PETSc C-based AIJ data
 * structures for parallel, sparse matrices. All overridden virtual
 * functions are documented in sparse_matrix.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief SparseMatrix interface to PETSc Mat.
 */
template <typename T>
class PetscMatrix final : public PetscMatrixBase<T>
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
  explicit
  PetscMatrix (const Parallel::Communicator & comm_in);

  /**
   * Constructor.  Creates a PetscMatrix assuming you already have a
   * valid Mat object.  In this case, m is NOT destroyed by the
   * PetscMatrix destructor when this object goes out of scope.  This
   * allows ownership of m to remain with the original creator, and to
   * simply provide additional functionality with the PetscMatrix.
   */
  explicit
  PetscMatrix (Mat m,
               const Parallel::Communicator & comm_in);

  /**
   * Constructor. Creates and initializes a PetscMatrix with the given
   * structure. See \p init(...) for a description of the parameters.
   */
  explicit
  PetscMatrix (const Parallel::Communicator & comm_in,
               const numeric_index_type m,
               const numeric_index_type n,
               const numeric_index_type m_l,
               const numeric_index_type n_l,
               const numeric_index_type n_nz=30,
               const numeric_index_type n_oz=10,
               const numeric_index_type blocksize=1);

  /**
   * This class manages a C-style struct (Mat) manually, so we
   * don't want to allow any automatic copy/move functions to be
   * generated, and we can't default the destructor.
   */
  PetscMatrix (PetscMatrix &&) = delete;
  PetscMatrix (const PetscMatrix &) = delete;
  PetscMatrix & operator= (PetscMatrix &&) = delete;
  virtual ~PetscMatrix ();

  PetscMatrix & operator= (const PetscMatrix &);
  virtual SparseMatrix<T> & operator= (const SparseMatrix<T> & v) override;

  /**
   * Initialize a PETSc matrix.
   *
   * \param m The global number of rows.
   * \param n The global number of columns.
   * \param m_l The local number of rows.
   * \param n_l The local number of columns.
   * \param n_nz The number of nonzeros in each row of the DIAGONAL portion of the local submatrix.
   * \param n_oz The number of nonzeros in each row of the OFF-DIAGONAL portion of the local submatrix.
   * \param blocksize Optional value indicating dense coupled blocks for systems with multiple variables all of the same type.
   */
  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type n_nz=30,
                     const numeric_index_type n_oz=10,
                     const numeric_index_type blocksize=1) override;

  /**
   * Initialize a PETSc matrix.
   *
   * \param m The global number of rows.
   * \param n The global number of columns.
   * \param m_l The local number of rows.
   * \param n_l The local number of columns.
   * \param n_nz array containing the number of nonzeros in each row of the DIAGONAL portion of the local submatrix.
   * \param n_oz Array containing the number of nonzeros in each row of the OFF-DIAGONAL portion of the local submatrix.
   * \param blocksize Optional value indicating dense coupled blocks for systems with multiple variables all of the same type.
   */
  void init (const numeric_index_type m,
             const numeric_index_type n,
             const numeric_index_type m_l,
             const numeric_index_type n_l,
             const std::vector<numeric_index_type> & n_nz,
             const std::vector<numeric_index_type> & n_oz,
             const numeric_index_type blocksize=1);

  virtual void init (ParallelType = PARALLEL) override;

  /**
   * Update the sparsity pattern based on \p dof_map, and set the matrix
   * to zero. This is useful in cases where the sparsity pattern changes
   * during a computation.
   */
  void update_preallocation_and_zero();

  /**
   * Reset matrix to use the original nonzero pattern provided by users
   */
  void reset_preallocation();

  virtual void zero () override;

  virtual std::unique_ptr<SparseMatrix<T>> zero_clone () const override;

  virtual std::unique_ptr<SparseMatrix<T>> clone () const override;

  virtual void zero_rows (std::vector<numeric_index_type> & rows, T diag_value = 0.0) override;

  virtual void flush () override;

  /**
   * Get the number of rows and columns owned by this process
   * @param i Row size
   * @param j Column size
   */
  void get_local_size (numeric_index_type & m, numeric_index_type & n) const;

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

  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & brows,
                                 const std::vector<numeric_index_type> & bcols) override;

  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & dof_indices) override
  { this->add_block_matrix (dm, dof_indices, dof_indices); }

  /**
   * Compute A += a*X for scalar \p a, matrix \p X.
   *
   * \note The matrices A and X need to have the same nonzero pattern,
   * otherwise \p PETSc will crash!
   *
   * \note It is advisable to not only allocate appropriate memory
   * with \p init(), but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.
   *
   * \note \p X will be closed, if not already done, before performing
   * any work.
   */
  virtual void add (const T a, const SparseMatrix<T> & X) override;

  /**
   * Compute Y = A*X for matrix \p X.
   * Set reuse = true if this->_mat and X have the same nonzero pattern as before,
   * and Y is obtained from a previous call to this function with reuse = false
   */
  virtual void matrix_matrix_mult (SparseMatrix<T> & X, SparseMatrix<T> & Y, bool reuse = false) override;

  /**
    * Add \p scalar* \p spm to the rows and cols of this matrix (A):
    * A(rows[i], cols[j]) += scalar * spm(i,j)
    */
  virtual
  void add_sparse_matrix (const SparseMatrix<T> & spm,
                          const std::map<numeric_index_type,numeric_index_type> & row_ltog,
                          const std::map<numeric_index_type,numeric_index_type> & col_ltog,
                          const T scalar) override;

  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const override;

  virtual Real l1_norm () const override;

  virtual Real linfty_norm () const override;

  /**
   * Print the contents of the matrix to the screen with the PETSc
   * viewer. This function only allows printing to standard out since
   * we have limited ourselves to one PETSc implementation for
   * writing.
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const override;

  virtual void print_matlab(const std::string & name = "") const override;

  virtual void print_petsc_binary(const std::string & filename) override;

  virtual void print_petsc_hdf5(const std::string & filename) override;

  virtual void read_petsc_binary(const std::string & filename) override;

  virtual void read_petsc_hdf5(const std::string & filename) override;

  virtual void get_diagonal (NumericVector<T> & dest) const override;

  virtual void get_transpose (SparseMatrix<T> & dest) const override;

  virtual
  void get_row(numeric_index_type i,
               std::vector<numeric_index_type> & indices,
               std::vector<T> & values) const override;

  /**
   * Similar to the `create_submatrix` function, this function creates a \p
   * submatrix which is defined by the indices given in the \p rows
   * and \p cols vectors.
   * Note: Both \p rows and \p cols can be unsorted;
   *       Use the above function for better efficiency if your indices are sorted;
   *       \p rows and \p cols can contain indices that are owned by other processors.
   */
   virtual void create_submatrix_nosort(SparseMatrix<T> & submatrix,
                                        const std::vector<numeric_index_type> & rows,
                                        const std::vector<numeric_index_type> & cols) const override;

  virtual void scale(const T scale) override;

protected:

  /**
   * This function either creates or re-initializes a matrix called \p
   * submatrix which is defined by the indices given in the \p rows
   * and \p cols vectors.
   *
   * This function is implemented in terms of MatGetSubMatrix().  The
   * \p reuse_submatrix parameter determines whether or not PETSc will
   * treat \p submatrix as one which has already been used (had memory
   * allocated) or as a new matrix.
   */
  virtual void _get_submatrix(SparseMatrix<T> & submatrix,
                              const std::vector<numeric_index_type> & rows,
                              const std::vector<numeric_index_type> & cols,
                              const bool reuse_submatrix) const override;

  void _petsc_viewer(const std::string & filename,
                     PetscViewerType viewertype,
                     PetscFileMode filemode);

  PetscMatrixType _mat_type;

private:
#ifdef LIBMESH_HAVE_CXX11_THREAD
  mutable std::mutex _petsc_matrix_mutex;
#else
  mutable Threads::spin_mutex _petsc_matrix_mutex;
#endif
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_MATRIX_H
