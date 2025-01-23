// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PETSC_MATRIX_BASE_H
#define LIBMESH_PETSC_MATRIX_BASE_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/wrapped_petsc.h"

// Petsc include files.
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscmat.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMatrix methods
#undef semiparallel_only
#undef exceptionless_semiparallel_only
#ifndef NDEBUG
#include <cstring>

#define semiparallel_only() do { if (this->initialized()) { const char * mytype; \
      LibmeshPetscCall(MatGetType(this->_mat,&mytype));                 \
      if (!strcmp(mytype, MATSEQAIJ))                                   \
        parallel_object_only(); } } while (0)
#define exceptionless_semiparallel_only() do { if (this->initialized()) { const char * mytype; \
      auto semiparallel_only_ierr = MatGetType(this->_mat,&mytype);     \
      libmesh_ignore(semiparallel_only_ierr);                           \
      if (!strcmp(mytype, MATSEQAIJ))                                   \
        exceptionless_parallel_object_only(); } } while (0)
#else
#define semiparallel_only()
#define exceptionless_semiparallel_only()
#endif


namespace libMesh
{

/**
 * This class provides a nice interface to the PETSc C-based data
 * structures for parallel, sparse matrices. All overridden virtual
 * functions are documented in sparse_matrix.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief SparseMatrix interface to PETSc Mat.
 */
template <typename T>
class PetscMatrixBase : public SparseMatrix<T>
{
public:
  explicit
  PetscMatrixBase (const Parallel::Communicator & comm_in);

  /**
   * Constructor.  Creates a PetscMatrixBase assuming you already have a
   * valid Mat object.  In this case, m may not be destroyed by the
   * PetscMatrixBase destructor when this object goes out of scope.  This
   * allows ownership of m to remain with the original creator, and to
   * simply provide additional functionality with the PetscMatrixBase.
   */
  explicit
  PetscMatrixBase (Mat m,
                   const Parallel::Communicator & comm_in,
                   bool destroy_on_exit = false);

  /**
   * This class manages a C-style struct (Mat) manually, so we
   * don't want to allow any automatic copy/move functions to be
   * generated, and we can't default the destructor.
   */
  PetscMatrixBase (PetscMatrixBase &&) = delete;
  PetscMatrixBase (const PetscMatrixBase &) = delete;
  PetscMatrixBase & operator= (PetscMatrixBase &&) = delete;
  virtual ~PetscMatrixBase ();

  virtual SolverPackage solver_package() override
  {
    return PETSC_SOLVERS;
  }

  /**
   * \returns The raw PETSc matrix pointer.
   *
   * \note This is generally not required in user-level code.
   *
   * \note Don't do anything crazy like calling MatDestroy() on
   * it, or very bad things will likely happen!
   */
  Mat mat () { libmesh_assert (_mat); return _mat; }

  /**
   * clear() is called from the destructor, so it should not throw.
   */
  virtual void clear() noexcept override;

  /**
   * If set to false, we don't delete the Mat on destruction and allow
   * instead for \p PETSc to manage it.
   */
  void set_destroy_mat_on_exit(bool destroy = true);

  /**
   * Swaps the internal data pointers of two PetscMatrices, no actual
   * values are swapped.
   */
  void swap (PetscMatrixBase<T> &);

  using SparseMatrix<T>::operator=;

  /**
   * Set the context (ourself) for \p _mat
   */
  void set_context();

  /**
   * @returns The context for \p mat if it exists, else a \p nullptr
   */
  static PetscMatrixBase<T> * get_context(Mat mat, const TIMPI::Communicator & comm);

  virtual numeric_index_type m () const override;

  virtual numeric_index_type local_m () const final;

  virtual numeric_index_type n () const override;

  /**
   * Get the number of columns owned by this process
   */
  virtual numeric_index_type local_n () const final;

  virtual numeric_index_type row_start () const override;

  virtual numeric_index_type row_stop () const override;

  virtual numeric_index_type col_start () const override;

  virtual numeric_index_type col_stop () const override;

  virtual void close () override;

  virtual bool closed() const override;

protected:

  /**
   * PETSc matrix datatype to store values.
   */
  Mat _mat;

  /**
   * This boolean value should only be set to \p false for the
   * constructor which takes a PETSc Mat object.
   */
  bool _destroy_mat_on_exit;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_MATRIX_BASE_H
