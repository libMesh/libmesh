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



#ifndef LIBMESH_MAPPED_SUBDOMAIN_PARTITIONER_H
#define LIBMESH_MAPPED_SUBDOMAIN_PARTITIONER_H

// Local Includes
#include "libmesh/partitioner.h"

namespace libMesh
{

/**
 * The \p MappedSubdomainPartitioner partitions the elements based on their
 * subdomain ids. The user must set up the values in the public
 * subdomain_to_proc map in order to specify which subdomains should
 * be partitioned onto which processors.  More than one subdomain can
 * be partitioned onto a given processor, but every subdomain must be
 * assigned to exactly 1 processor.
 *
 * \author John W. Peterson
 * \date 2017
 * \brief Partitions elements based on user-defined mapping from subdomain ids -> processor ids.
 */
class MappedSubdomainPartitioner : public Partitioner
{
public:

  /**
   * Constructor.
   */
  MappedSubdomainPartitioner () {}

  /**
   * Creates a new partitioner of this type and returns it in
   * a \p UniquePtr.
   */
  virtual UniquePtr<Partitioner> clone () const libmesh_override
  {
    return UniquePtr<Partitioner>(new MappedSubdomainPartitioner());
  }

  /**
   * Before calling partition() or partition_range(), the user must
   * assign all the Mesh subdomains to certain processors by adding
   * them to this std::map.  For example:
   * subdomain_to_proc[1] = 0;
   * subdomain_to_proc[2] = 0;
   * subdomain_to_proc[3] = 0;
   * subdomain_to_proc[4] = 1;
   * subdomain_to_proc[5] = 1;
   * subdomain_to_proc[6] = 2;
   * Would partition the mesh onto three processors, with subdomains
   * 1, 2, and 3 on processor 0, subdomains 4 and 5 on processor 1,
   * and subdomain 6 on processor 2.
   */
  std::map<subdomain_id_type, processor_id_type> subdomain_to_proc;

  /**
   * Called by the SubdomainPartitioner to partition elements in the range (it, end).
   */
  virtual void partition_range(MeshBase & mesh,
                               MeshBase::element_iterator it,
                               MeshBase::element_iterator end,
                               const unsigned int n) libmesh_override;

protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) libmesh_override;
};

} // namespace libMesh

#endif  // LIBMESH_MAPPED_SUBDOMAIN_PARTITIONER_H
