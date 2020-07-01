// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MATRIX_H
#define LIBMESH_MATRIX_H


// Local includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"

namespace libMesh
{
/**
 * Generic matrix. Currently is nothing more than a stand-in allowing for covariant return types
 * from things like \p ImplicitSystem and \p EigenSystem
 *
 * \author Alexander Lindsay
 * \date 2020
 */
class Matrix : public ReferenceCountedObject<Matrix>,
               public ParallelObject
{
public:
  /**
   * Constructor; does nothing.
   */
  Matrix (const Parallel::Communicator & comm_in) : ParallelObject(comm_in) {}

  /**
   * Destructor.
   */
  virtual ~Matrix () = default;
};
} // namespace libMesh


#endif // LIBMESH_MATRIX_H
