// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/parallel_ghost_sync.h"

namespace libMesh
{

SyncNodalPositions::SyncNodalPositions(MeshBase & m)
  : mesh(m)
{}



void SyncNodalPositions::gather_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data) const
{
  data.resize(ids.size());

  // Gather (x,y,z) data for all node IDs in the ids vector
  for (std::size_t i=0; i<ids.size(); ++i)
    {
      // Look for this node in the mesh
      const Point & pt = mesh.point(ids[i]);

      // Store this node's position in the data array.
      // This should call Point::op=
      data[i] = pt;
    } // end for
} // gather_data()



void SyncNodalPositions::act_on_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data) const
{
  for (std::size_t i=0; i<ids.size(); ++i)
    {

      // Get a pointer to the node whose position is to be updated.
      Node & node = mesh.node_ref(ids[i]);

      // Update this node's position.  Should call Point::op=
      node = data[i];
    } // end for
} // act_on_data()

} // namespace
