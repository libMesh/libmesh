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

#ifndef LIBMESH_PETSC_MATRIX_SHELL_MATRIX_H
#define LIBMESH_PETSC_MATRIX_SHELL_MATRIX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_shell_matrix.h"

namespace libMesh
{

/**
 * This class allows to use a PETSc shell matrix as a PetscMatrix
 *
 * \author Alexander Lindsay
 * \date 2024
 */
template <typename T>
class PetscMatrixShellMatrix : public PetscMatrixBase<T>
{
public:
  explicit PetscMatrixShellMatrix(const Parallel::Communicator & comm_in);

  virtual void init(const numeric_index_type m,
                    const numeric_index_type n,
                    const numeric_index_type m_l,
                    const numeric_index_type n_l,
                    const numeric_index_type = 30,
                    const numeric_index_type = 10,
                    const numeric_index_type blocksize = 1) override;

  /**
   * Initialize this matrix using the sparsity structure computed by \p dof_map.
   * @param type The serial/parallel/ghosted type of the matrix
   */
  virtual void init(ParallelType = PARALLEL) override;

  virtual SparseMatrix<T> & operator=(const SparseMatrix<T> &) override;

private:
  // Make this private because we mark as initialized after we've done our initialization, and we
  // don't want derived classes to mistakenly register their data as initialized (or not)
  using PetscMatrixBase<T>::_is_initialized;

  /// Whether to omit constrained degrees of freedom
  const bool _omit_constrained_dofs;

  friend void init_shell_mat<PetscMatrixShellMatrix<T>>(PetscMatrixShellMatrix<T> & obj);
  friend void init_shell_mat<PetscMatrixShellMatrix<T>>(PetscMatrixShellMatrix<T> & obj,
                                                        const numeric_index_type m,
                                                        const numeric_index_type n,
                                                        const numeric_index_type m_l,
                                                        const numeric_index_type n_l,
                                                        const numeric_index_type blocksize_in);
};

//-----------------------------------------------------------------------
// PetscMatrixShellMatrix inline members
template <typename T>
PetscMatrixShellMatrix<T>::PetscMatrixShellMatrix(const Parallel::Communicator & comm_in)
  : PetscMatrixBase<T>(comm_in), _omit_constrained_dofs(false)
{
}

template <typename T>
SparseMatrix<T> &
PetscMatrixShellMatrix<T>::operator=(const SparseMatrix<T> &)
{
  libmesh_error();
}

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_SPARSE_SHELL_MATRIX_H
