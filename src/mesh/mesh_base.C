// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute and/or
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

#include "libmesh/mesh_base_impl.h"


namespace libMesh
{
template <>
void MeshBaseTempl<Real>::partition (const unsigned int n_parts)
{
  // If we get here and we have unpartitioned elements, we need that
  // fixed.
  if (this->n_unpartitioned_elem() > 0)
    {
      libmesh_assert (partitioner().get());
      libmesh_assert (this->is_serial());
      partitioner()->partition (*this, n_parts);
    }
  // A nullptr partitioner or a skip_partitioning(true) call or a
  // skip_noncritical_partitioning(true) call means don't repartition;
  // skip_noncritical_partitioning() checks all these.
  else if (!skip_noncritical_partitioning())
    {
      partitioner()->partition (*this, n_parts);
    }
  else
    {
      // Adaptive coarsening may have "orphaned" nodes on processors
      // whose elements no longer share them.  We need to check for
      // and possibly fix that.
      MeshTools::correct_node_proc_ids(*this);

      // Make sure locally cached partition count is correct
      this->recalculate_n_partitions();

      // Make sure any other locally cached data is correct
      this->update_post_partitioning();
    }
}


template class MeshBaseTempl<Real>;
}
