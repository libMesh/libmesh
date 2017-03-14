// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SUBDOMAIN_PARTITIONER_H
#define LIBMESH_SUBDOMAIN_PARTITIONER_H

// Local Includes
#include "libmesh/partitioner.h"

namespace libMesh
{

/**
 * The \p SubdomainPartitioner partitions the elements in "chunks" of
 * user-specified subdomain ids.  Once all the chunks are partitioned,
 * the overall mesh partitioning is simply the union of the chunk
 * partitionings. For example, if the "chunks" vector is given by:
 * chunks[0] = {1, 2, 4}
 * chunks[1] = {3, 7, 8}
 * chunks[2] = {5, 6}
 * then we will call the internal Partitioner three times, once for
 * subodmains 1, 2, and 4, once for subdomains 3, 7, and 8, and once
 * for subdomains 5 and 6.
 *
 * Note: this Partitioner may produce highly non-optimal communication
 * patterns and is likely to place geometrically disjoint sets of
 * elements on the same processor. Its intended use is to help
 * facilitate load balancing. That is, if the user knows that certain
 * subdomains (or groups of subdomains) are more expensive to compute
 * than others, he/she can ensure that they are partitioned more or
 * less evenly among the available processors by specifying them
 * together in a single entry of the "chunk" vector.
 *
 * \author John W. Peterson
 * \date 2017
 * \brief Independently partitions chunks of subdomains and combines the results.
 */
class SubdomainPartitioner : public Partitioner
{
public:

  /**
   * Constructor, default initializes the internal Partitioner object
   * to a MetisPartitioner so the class is usable, although this type
   * can be customized later.
   */
  SubdomainPartitioner ();

  /**
   * Creates a new partitioner of this type and returns it in
   * a \p UniquePtr.
   */
  virtual UniquePtr<Partitioner> clone () const libmesh_override
  {
    return UniquePtr<Partitioner>(new SubdomainPartitioner());
  }

  /**
   * Each entry of "chunks" represents a set of subdomains which are
   * to be partitioned together.  The internal Partitioner will be
   * called once for each entry of chunks, and the resulting
   * partitioning will simply be the union of all these partitionings.
   */
  std::vector<std::set<subdomain_id_type> > chunks;

  /**
   * Get a reference to the Partitioner used internally by the
   * SubdomainPartitioner.  Note: the internal Partitioner cannot also
   * be a SubdomainPartitioner, otherwise an infinite loop will
   * result.  To have this class use e.g. the space-filling curve
   * Partitioner internally, one could do:
   *
   * SubdomainPartitioner sp;
   * sp.internal_partitioner().reset(new SFCPartitioner);
   */
  UniquePtr<Partitioner> & internal_partitioner() { return _internal_partitioner; }

protected:
  /**
   * The internal Partitioner we use. Public access via the
   * internal_partitioner() member function.
   */
  UniquePtr<Partitioner> _internal_partitioner;

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) libmesh_override;
};

} // namespace libMesh

#endif  // LIBMESH_SUBDOMAIN_PARTITIONER_H
