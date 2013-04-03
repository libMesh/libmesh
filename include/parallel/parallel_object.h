// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/parallel.h"

// Macro to identify and debug functions which should only be called in
// parallel on every processor at once
#undef parallel_object_only
#ifndef NDEBUG
  #define parallel_object_only() do { \
    libmesh_assert(this->communicator().verify(std::string(__FILE__).size())); \
    libmesh_assert(this->communicator().verify(std::string(__FILE__))); \
    libmesh_assert(this->communicator().verify(__LINE__)); } while (0)
#else
  #define parallel_object_only()  ((void) 0)
#endif


namespace libMesh
{
  /**
   * This class forms the base class for all other classes
   * that are expected to be implemented in paralel. Each
   * \p ParalelObject *requires* a \p Parallel::Communicator object
   * for construction.
   *
   * \author Benjamin S. Kirk, 2013.
   */
  class ParallelObject
  {
  public:

    /**
     * Constructor. Requires a reference to the communicator
     * that defines the object's parallel decomposition.
     */
    ParallelObject (const Parallel::Communicator &comm) :
      _communicator(comm)
    {}

    /**
     * Copy Constructor.
     */
    ParallelObject (const ParallelObject &other) :
      _communicator(other._communicator)
    {}

    /**
     * Destructor.  Virtual because we are a base class.
     */
    virtual ~ParallelObject () {};

    /**
     * @returns a reference to the \p Parallel::Communicator object
     * used by this mesh.
     */
    const Parallel::Communicator & communicator () const
    { return _communicator; }

    /**
     * @returns the number of processors in the group.
     */
    processor_id_type n_processors () const
    { return libmesh_cast_int<processor_id_type>(_communicator.size()); }

    /**
     * @returns the rank of this processor in the group.
     */
    processor_id_type processor_id () const
    { return libmesh_cast_int<processor_id_type>(_communicator.rank()); }


  protected:

    const Parallel::Communicator &_communicator;
  };
} // namespace libMesh

#endif // LIBMESH_PARALLEL_OBJECT_H
