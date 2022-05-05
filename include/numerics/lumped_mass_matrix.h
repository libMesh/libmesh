// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_LUMPED_MASS_MATRIX_H
#define LIBMESH_LUMPED_MASS_MATRIX_H

#include "libmesh/diagonal_matrix.h"

namespace libMesh
{

/**
 * Template class used to construct a lumped mass matrix. Potentially also useful for computing
 * quantities relevant to overall system scaling. Any time the add method is called on this class we
 * sum into the row index \p i with the absolute value of the provided value
 *
 * \author Alexander D. Lindsay
 * \date 2021
 */
template <typename T>
class LumpedMassMatrix : public DiagonalMatrix<T>
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
   * You have to initialize the matrix before usage with
   * \p init(...).
   */
  explicit LumpedMassMatrix(const Parallel::Communicator & comm);

  /**
   * unique pointers can be moved but not copied
   */
  LumpedMassMatrix(LumpedMassMatrix &&) = default;
  LumpedMassMatrix & operator=(LumpedMassMatrix &&) = default;

  virtual std::unique_ptr<SparseMatrix<T>> zero_clone() const override;

  virtual std::unique_ptr<SparseMatrix<T>> clone() const override;

  virtual void set(const numeric_index_type i,
                   const numeric_index_type j,
                   const T value) override;

  virtual void add(const numeric_index_type i,
                   const numeric_index_type j,
                   const T value) override;
  virtual void add(const T a, const SparseMatrix<T> & X) override;

protected:
  /**
   * Copy contents from vec into underlying diagonal storage
   */
  LumpedMassMatrix & operator=(const NumericVector<T> & vec);

  /**
   * Move contents from vec into underlying diagonal storage
   */
  LumpedMassMatrix & operator=(NumericVector<T> && vec);
};

} // namespace libMesh

#endif // LIBMESH_LUMPED_MASS_MATRIX_H
