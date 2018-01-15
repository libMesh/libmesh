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



#ifndef LIBMESH_TENSOR_SHELL_MATRIX_H
#define LIBMESH_TENSOR_SHELL_MATRIX_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{

/**
 * Shell matrix that is given by a tensor product of two vectors,
 * i.e. A = v*w^T. All overridden virtual functions are documented in
 * shell_matrix.h.
 *
 * \author Tim Kroeger
 * \date 2008
 */
template <typename T>
class TensorShellMatrix : public ShellMatrix<T>
{
public:
  /**
   * Constructor; takes references to the two vectors as arguments.
   * The vectors themselves have to be stored elsewhere.
   */
  TensorShellMatrix (const NumericVector<T> & v,
                     const NumericVector<T> & w);

  /**
   * Destructor.
   */
  virtual ~TensorShellMatrix ();

  virtual numeric_index_type m () const libmesh_override;

  virtual numeric_index_type n () const libmesh_override;

  virtual void vector_mult (NumericVector<T> & dest,
                            const NumericVector<T> & arg) const libmesh_override;

  virtual void vector_mult_add (NumericVector<T> & dest,
                                const NumericVector<T> & arg) const libmesh_override;

  virtual void get_diagonal (NumericVector<T> & dest) const libmesh_override;

protected:
  /**
   * The column vector.
   */
  const NumericVector<T> & _v;

  /**
   * The row vector.
   */
  const NumericVector<T> & _w;
};



//-----------------------------------------------------------------------
// TensorShellMatrix inline members
template <typename T>
inline
TensorShellMatrix<T>::TensorShellMatrix (const NumericVector<T> & v,
                                         const NumericVector<T> & w):
  ShellMatrix<T>(v.comm()),
  _v(v),
  _w(w)
{}



template <typename T>
inline
TensorShellMatrix<T>::~TensorShellMatrix ()
{}



template <typename T>
inline
numeric_index_type TensorShellMatrix<T>::m () const
{
  return _v.size();
}



template <typename T>
inline
numeric_index_type TensorShellMatrix<T>::n () const
{
  return _w.size();
}


} // namespace libMesh


#endif // LIBMESH_TENSOR_SHELL_MATRIX_H
