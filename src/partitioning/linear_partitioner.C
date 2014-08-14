// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh/mesh_base.h"
#include "libmesh/linear_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"

namespace libMesh
{


// ------------------------------------------------------------
// LinearPartitioner implementation
void LinearPartitioner::_do_partition (MeshBase& mesh,
                                       const unsigned int n)
{
  libmesh_assert_greater (n, 0);

  // Check for an easy return
  if (n == 1)
    {
      this->single_partition (mesh);
      return;
    }

  // Create a simple linear partitioning
  {
    START_LOG ("partition()", "LinearPartitioner");

    const dof_id_type n_active_elem = mesh.n_active_elem();
    const dof_id_type blksize       = n_active_elem/n;

    dof_id_type e = 0;

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      {
        if ((e/blksize) < n)
          {
            Elem *elem = *elem_it;
            elem->processor_id() =
              cast_int<processor_id_type>(e/blksize);
          }
        else
          {
            Elem *elem = *elem_it;
            elem->processor_id() = 0;
            elem = elem->parent();
          }

        e++;
      }

    STOP_LOG ("partition()", "LinearPartitioner");
  }
}

} // namespace libMesh
