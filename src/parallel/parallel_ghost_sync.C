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



#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/int_range.h"

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
  for (auto i : index_range(ids))
    {
      // Look for this node in the mesh
      const Point & pt = mesh.point(ids[i]);

      // Store this node's position in the data array.
      // This should call Point::op=
      data[i] = pt;
    } // end for
} // gather_data()



void SyncNodalPositions::act_on_data (const std::vector<dof_id_type> & ids,
                                      const std::vector<datum> & data) const
{
  for (auto i : index_range(ids))
    {
      // Get a pointer to the node whose position is to be updated.
      Node & node = mesh.node_ref(ids[i]);

      // Update this node's position.  Should call Point::op=
      node = data[i];
    } // end for
} // act_on_data()


SyncSubdomainIds::SyncSubdomainIds(MeshBase & m)
  : mesh(m)
{}

void SyncSubdomainIds::gather_data (const std::vector<dof_id_type> & ids,
                                        std::vector<datum> & ids_out) const
{
  ids_out.reserve(ids.size());

  for (const auto & id : ids)
    {
      Elem & elem = mesh.elem_ref(id);
      ids_out.push_back(elem.subdomain_id());
    }
}

void SyncSubdomainIds::act_on_data (const std::vector<dof_id_type> & ids,
                              const std::vector<datum> & subdomain_ids) const
{
  for (auto i : index_range(ids))
    {
      Elem & elem = mesh.elem_ref(ids[i]);
      elem.subdomain_id()=subdomain_ids[i];
    }
}

} // namespace
