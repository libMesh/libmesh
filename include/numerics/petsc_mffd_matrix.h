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

#ifndef LIBMESH_PETSC_MFFD_MATRIX_H
#define LIBMESH_PETSC_MFFD_MATRIX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_matrix_base.h"

namespace libMesh
{

/**
 * This class allows to use a PETSc shell matrix as a PetscMatrix
 *
 * \author Alexander Lindsay
 * \date 2024
 */
template <typename T>
class PetscMFFDMatrix : public PetscMatrixBase<T>
{
public:
  /**
   * Constructor.  Creates a PetscMFFDMatrix assuming you already have a
   * valid Mat object.  In this case, m is NOT destroyed by the
   * PetscMFFDMatrix destructor when this object goes out of scope.  This
   * allows ownership of m to remain with the original creator, and to
   * simply provide additional functionality with the PetscMFFDMatrix.
   */
  explicit PetscMFFDMatrix(Mat m, const Parallel::Communicator & comm_in);

  explicit PetscMFFDMatrix(const Parallel::Communicator & comm_in);

  PetscMFFDMatrix & operator=(Mat m);

  virtual void init(const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type = 30,
                    const numeric_index_type = 10,
                    const numeric_index_type = 1) override;

  /**
   * Initialize this matrix using the sparsity structure computed by \p dof_map.
   * @param type The serial/parallel/ghosted type of the matrix
   */
  virtual void init(ParallelType = PARALLEL) override;

  virtual SparseMatrix<T> & operator=(const SparseMatrix<T> &) override;

  virtual void zero() override;
  virtual std::unique_ptr<SparseMatrix<T>> zero_clone() const override;
  virtual std::unique_ptr<SparseMatrix<T>> clone() const override;
  virtual void set(const numeric_index_type i, const numeric_index_type j, const T value) override;
  virtual void add(const numeric_index_type i, const numeric_index_type j, const T value) override;
  virtual void add_matrix(const DenseMatrix<T> & dm,
                          const std::vector<numeric_index_type> & rows,
                          const std::vector<numeric_index_type> & cols) override;
  virtual void add_matrix(const DenseMatrix<T> & dm,
                          const std::vector<numeric_index_type> & dof_indices) override;
  virtual void add(const T a, const SparseMatrix<T> & X) override;
  virtual T operator()(const numeric_index_type i, const numeric_index_type j) const override;
  virtual Real l1_norm() const override;
  virtual Real linfty_norm() const override;
  virtual void print_personal(std::ostream & os = libMesh::out) const override;
  virtual void get_diagonal(NumericVector<T> & dest) const override;
  virtual void get_transpose(SparseMatrix<T> & dest) const override;
  virtual void get_row(numeric_index_type i,
                       std::vector<numeric_index_type> & indices,
                       std::vector<T> & values) const override;
};

template <typename T>
PetscMFFDMatrix<T>::PetscMFFDMatrix(Mat m, const Parallel::Communicator & comm_in)
  : PetscMatrixBase<T>(m, comm_in)
{
}

template <typename T>
PetscMFFDMatrix<T>::PetscMFFDMatrix(const Parallel::Communicator & comm_in)
  : PetscMatrixBase<T>(comm_in)
{
}

template <typename T>
PetscMFFDMatrix<T> &
PetscMFFDMatrix<T>::operator=(Mat m)
{
  this->_mat = m;
  return *this;
}

template <typename T>
void
PetscMFFDMatrix<T>::init(const numeric_index_type,
                         const numeric_index_type,
                         const numeric_index_type,
                         const numeric_index_type,
                         const numeric_index_type,
                         const numeric_index_type,
                         const numeric_index_type)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::init(ParallelType)
{
  libmesh_error();
}

template <typename T>
SparseMatrix<T> &
PetscMFFDMatrix<T>::operator=(const SparseMatrix<T> &)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::zero()
{
  libmesh_error();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
PetscMFFDMatrix<T>::zero_clone() const
{
  libmesh_error();
}

template <typename T>
std::unique_ptr<SparseMatrix<T>>
PetscMFFDMatrix<T>::clone() const
{
  libmesh_not_implemented();
}

template <typename T>
void
PetscMFFDMatrix<T>::set(const numeric_index_type, const numeric_index_type, const T)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::add(const numeric_index_type, const numeric_index_type, const T)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::add_matrix(const DenseMatrix<T> &,
                               const std::vector<numeric_index_type> &,
                               const std::vector<numeric_index_type> &)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::add_matrix(const DenseMatrix<T> &, const std::vector<numeric_index_type> &)
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::add(const T, const SparseMatrix<T> &)
{
  libmesh_error();
}

template <typename T>
T
PetscMFFDMatrix<T>::operator()(const numeric_index_type, const numeric_index_type) const
{
  libmesh_error();
}

template <typename T>
Real
PetscMFFDMatrix<T>::l1_norm() const
{
  libmesh_error();
}

template <typename T>
Real
PetscMFFDMatrix<T>::linfty_norm() const
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::print_personal(std::ostream &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::get_diagonal(NumericVector<T> &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::get_transpose(SparseMatrix<T> &) const
{
  libmesh_error();
}

template <typename T>
void
PetscMFFDMatrix<T>::get_row(numeric_index_type,
                            std::vector<numeric_index_type> &,
                            std::vector<T> &) const
{
  libmesh_error();
}

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_SPARSE_SHELL_MATRIX_H
