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



#ifndef LIBMESH_PARALLEL_OBJECT_H
#define LIBMESH_PARALLEL_OBJECT_H


// libMesh includes
#include "libmesh/libmesh_common.h" // libmesh_dbg_var
#include "libmesh/parallel_only.h"

// TIMPI includes
#include "timpi/communicator.h"


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once
#undef parallel_object_only
#undef exceptionless_parallel_object_only
#ifndef NDEBUG
#define parallel_object_only() libmesh_parallel_only(this->comm())
#define exceptionless_parallel_object_only() libmesh_exceptionless_parallel_only(this->comm())
#else
#define parallel_object_only()  ((void) 0)
#define exceptionless_parallel_object_only()  ((void) 0)
#endif


namespace libMesh
{
/**
 * \brief An object whose state is distributed along a set of processors.
 *
 * This class forms the base class for all other classes
 * that are expected to be implemented in parallel. Each
 * \p ParallelObject *requires* a \p Parallel::Communicator object
 * for construction.
 *
 * \author Benjamin S. Kirk
 * \date 2013
 */
class ParallelObject
{
public:

  /**
   * Constructor. Requires a reference to the communicator
   * that defines the object's parallel decomposition.
   */
  ParallelObject (const Parallel::Communicator & comm_in) :
    _communicator(comm_in)
  {}

  /**
   * Copy Constructor.
   */
  ParallelObject (const ParallelObject & other) :
    _communicator(other._communicator)
  {}

  /**
   * "Assignment" operator.  Simply asserts our references
   * are identical because this is the only thing that makes
   * sense
   */
  ParallelObject & operator= (const ParallelObject & libmesh_dbg_var(other))
  {
    libmesh_assert_equal_to (&_communicator, &other._communicator);
    return *this;
  }

  /**
   * Destructor.  Virtual because we are a base class.
   */
  virtual ~ParallelObject () = default;

  /**
   * \returns A reference to the \p Parallel::Communicator object
   * used by this mesh.
   */
  const Parallel::Communicator & comm () const
  { return _communicator; }

  /**
   * \returns The number of processors in the group.
   */
  processor_id_type n_processors () const
  {
    processor_id_type returnval =
      cast_int<processor_id_type>(_communicator.size());
    libmesh_assert(returnval); // We never have an empty comm
    return returnval;
  }

  /**
   * \returns The rank of this processor in the group.
   */
  processor_id_type processor_id () const
  { return cast_int<processor_id_type>(_communicator.rank()); }


protected:

  const Parallel::Communicator & _communicator;
};
} // namespace libMesh

#endif // LIBMESH_PARALLEL_OBJECT_H
