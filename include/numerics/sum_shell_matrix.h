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



#ifndef LIBMESH_SUM_SHELL_MATRIX_H
#define LIBMESH_SUM_SHELL_MATRIX_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/shell_matrix.h"

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * This class combines any number of shell matrices to a single shell
 * matrix, acting as the sum of the matrices.
 *
 * \author Tim Kroeger
 * \date 2008
 */
template <typename T>
class SumShellMatrix : public ShellMatrix<T>
{
public:
  /**
   * Constructor; initializes an empty sum.  Note that an empty sum is
   * not a valid object in that a call to \p m() or \p n() will result
   * in an error.  However, an empty sum is allowed to be multiplied
   * with a vector and will give the expected result.
   */
  SumShellMatrix (const Parallel::Communicator & comm_in
                  LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Constructor that passes a vector of shell matrices.
   */
  explicit
  SumShellMatrix (const std::vector<ShellMatrix<T> *> & mat,
                  const Parallel::Communicator & comm_in
                  LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~SumShellMatrix ();

  /**
   * @returns \p m, the row-dimension of the matrix where the marix is
   * \f$ M \times N \f$.
   */
  virtual numeric_index_type m () const libmesh_override;

  /**
   * @returns \p n, the column-dimension of the matrix where the marix
   * is \f$ M \times N \f$.
   */
  virtual numeric_index_type n () const libmesh_override;

  /**
   * Multiplies the matrix with \p arg and stores the result in \p
   * dest.
   */
  virtual void vector_mult (NumericVector<T> & dest,
                            const NumericVector<T> & arg) const libmesh_override;

  /**
   * Multiplies the matrix with \p arg and adds the result to \p dest.
   */
  virtual void vector_mult_add (NumericVector<T> & dest,
                                const NumericVector<T> & arg) const libmesh_override;

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T> & dest) const libmesh_override;

  /**
   * A vector of the summands.
   */
  std::vector<ShellMatrix<T> *> matrices;
};



//-----------------------------------------------------------------------
// SumShellMatrix inline members
template <typename T>
inline
SumShellMatrix<T>::SumShellMatrix (const Parallel::Communicator & comm_in):
  ShellMatrix<T>(comm_in),
  matrices()
{}



template <typename T>
inline
SumShellMatrix<T>::SumShellMatrix (const std::vector<ShellMatrix<T> *> & mat,
                                   const Parallel::Communicator & comm_in):
  ShellMatrix<T>(comm_in),
  matrices(mat)
{}



template <typename T>
inline
SumShellMatrix<T>::~SumShellMatrix ()
{}


} // namespace libMesh


#endif // LIBMESH_SUM_SHELL_MATRIX_H
