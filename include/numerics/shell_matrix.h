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



#ifndef LIBMESH_SHELL_MATRIX_H
#define LIBMESH_SHELL_MATRIX_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"

namespace libMesh
{

// forward declarations
template <typename T> class NumericVector;


/**
 * Generic shell matrix, i.e. a matrix that does not define anything
 * but its action on a vector.  This class contains pure virtual
 * members that must be overloaded in derived classes.
 *
 * \author Tim Kroeger
 * \date 2008
 */
template <typename T>
class ShellMatrix : public ReferenceCountedObject<ShellMatrix<T> >,
                    public ParallelObject
{
public:
  /**
   * Constructor; does nothing.
   */
  ShellMatrix (const Parallel::Communicator & comm_in
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~ShellMatrix ();

  /**
   * @returns \p m, the row-dimension of the matrix where the marix is
   * \f$ M \times N \f$.
   */
  virtual numeric_index_type m () const = 0;

  /**
   * @returns \p n, the column-dimension of the matrix where the marix
   * is \f$ M \times N \f$.
   */
  virtual numeric_index_type n () const = 0;

  /**
   * Multiplies the matrix with \p arg and stores the result in \p
   * dest.
   */
  virtual void vector_mult (NumericVector<T> & dest,
                            const NumericVector<T> & arg) const = 0;

  /**
   * Multiplies the matrix with \p arg and adds the result to \p dest.
   */
  virtual void vector_mult_add (NumericVector<T> & dest,
                                const NumericVector<T> & arg) const = 0;

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T> & dest) const = 0;
};



//-----------------------------------------------------------------------
// ShellMatrix inline members
template <typename T>
inline
ShellMatrix<T>::ShellMatrix (const Parallel::Communicator & comm_in) :
  ParallelObject(comm_in)
{}



template <typename T>
inline
ShellMatrix<T>::~ShellMatrix ()
{}


} // namespace libMesh


#endif // LIBMESH_SHELL_MATRIX_H
