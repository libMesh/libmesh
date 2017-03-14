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
#include "libmesh/mesh_base.h"
#include "libmesh/mapped_subdomain_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"

namespace libMesh
{


void MappedSubdomainPartitioner::partition_range(MeshBase & /*mesh*/,
                                                 MeshBase::element_iterator it,
                                                 MeshBase::element_iterator end,
                                                 const unsigned int n)
{
  libmesh_assert_greater (n, 0);

  // Check for an easy return.  If the user has not specified any
  // entries in subdomain_to_proc, we just do a single partitioning.
  if ((n == 1) || subdomain_to_proc.empty())
    {
      this->single_partition_range (it, end);
      return;
    }

  // Now actually do the partitioning.
  LOG_SCOPE ("partition_range()", "MappedSubdomainPartitioner");

  for ( ; it != end; ++it)
    {
      Elem * elem = *it;

      subdomain_id_type sbd_id = elem->subdomain_id();

      // Find which processor id corresponds to this element's subdomain id.
      std::map<subdomain_id_type, processor_id_type>::iterator
        it = subdomain_to_proc.find(sbd_id);

      // If an element has a subdomain we don't know how to
      // partition, throw an error.
      if (it == subdomain_to_proc.end())
        libmesh_error_msg("Could not find processor id corresponding to subdomain id " << sbd_id);

      elem->processor_id() = it->second;
    }
}



void MappedSubdomainPartitioner::_do_partition (MeshBase & mesh,
                                                const unsigned int n)
{
  this->partition_range(mesh,
                        mesh.active_elements_begin(),
                        mesh.active_elements_end(),
                        n);
}

} // namespace libMesh
