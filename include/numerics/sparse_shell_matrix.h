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



#ifndef LIBMESH_SPARSE_SHELL_MATRIX_H
#define LIBMESH_SPARSE_SHELL_MATRIX_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

/**
 * This class allows to use any SparseMatrix object as a shell matrix.
 *
 * \author Tim Kroeger
 * \date 2008
 */
template <typename T>
class SparseShellMatrix : public ShellMatrix<T>
{
public:
  /**
   * Constructor; takes references to the sparse matrix.  The sparse
   * matrix itself has to be stored elsewhere.
   */
  explicit
  SparseShellMatrix (const SparseMatrix<T> & new_m);

  /**
   * Destructor.
   */
  virtual ~SparseShellMatrix ();

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

protected:
  /**
   * The sparse matrix.
   */
  const SparseMatrix<T> & _m;
};



//-----------------------------------------------------------------------
// SparseShellMatrix inline members
template <typename T>
inline
SparseShellMatrix<T>::SparseShellMatrix (const SparseMatrix<T> & new_m):
  ShellMatrix<T>(new_m.comm()),
  _m(new_m)
{}



template <typename T>
inline
SparseShellMatrix<T>::~SparseShellMatrix ()
{}



template <typename T>
inline
numeric_index_type SparseShellMatrix<T>::m () const
{
  return _m.m();
}



template <typename T>
inline
numeric_index_type SparseShellMatrix<T>::n () const
{
  return _m.n();
}



template <typename T>
inline
void SparseShellMatrix<T>::get_diagonal(NumericVector<T> & dest) const
{
  _m.get_diagonal(dest);
}


} // namespace libMesh


#endif // LIBMESH_SPARSE_SHELL_MATRIX_H
