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



// Local Includes
#include "libmesh/subdomain_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/metis_partitioner.h"

namespace libMesh
{

SubdomainPartitioner::SubdomainPartitioner () :
  _internal_partitioner(new MetisPartitioner)
{}


void SubdomainPartitioner::_do_partition (MeshBase & mesh,
                                          const unsigned int n)
{
  libmesh_assert_greater (n, 0);

  // Check for an easy return.  If the user has not specified any
  // entries in the chunks vector, we just do a single partitioning.
  if ((n == 1) || chunks.empty())
    {
      this->single_partition (mesh);
      return;
    }

  // Now actually do the partitioning.
  LOG_SCOPE ("_do_partition()", "SubdomainPartitioner");

  // For each chunk, construct an iterator range for the set of
  // subdomains in question, and pass it to the internal Partitioner.
  for (std::size_t c=0; c<chunks.size(); ++c)
    {
      MeshBase::element_iterator
        it = mesh.active_subdomain_set_elements_begin(chunks[c]),
        end = mesh.active_subdomain_set_elements_end(chunks[c]);

      _internal_partitioner->partition_range(mesh, it, end, n);
    }
}

} // namespace libMesh
