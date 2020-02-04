// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/topology_map_impl.h"

namespace libMesh
{


dof_id_type TopologyMap::find(const std::vector<std::pair<dof_id_type, dof_id_type>> & bracketing_nodes) const
{
  dof_id_type new_node_id = DofObject::invalid_id;

  for (auto pair : bracketing_nodes)
    {
      const dof_id_type lower_id = std::min(pair.first, pair.second);
      const dof_id_type upper_id = std::max(pair.first, pair.second);

      const dof_id_type possible_new_node_id =
        this->find(lower_id, upper_id);

      if (possible_new_node_id != DofObject::invalid_id)
        {
          // If we found a node already, but we're still here, it's to
          // debug map consistency: we'd better always find the same
          // node
          if (new_node_id != DofObject::invalid_id)
            libmesh_assert_equal_to (new_node_id, possible_new_node_id);

          new_node_id = possible_new_node_id;
        }

      // If we're not debugging map consistency then we can quit as
      // soon as we find a node
#ifdef NDEBUG
      if (new_node_id != DofObject::invalid_id)
        break;
#endif
    }
  return new_node_id;
}


dof_id_type TopologyMap::find(dof_id_type bracket_node1,
                              dof_id_type bracket_node2) const
{
  const dof_id_type lower_id = std::min(bracket_node1, bracket_node2);
  const dof_id_type upper_id = std::max(bracket_node1, bracket_node2);

  map_type::const_iterator it =
    _map.find(std::make_pair(lower_id, upper_id));

  if (it == _map.end())
    return DofObject::invalid_id;

  libmesh_assert_not_equal_to (it->second, DofObject::invalid_id);

  return it->second;
}

TOPO_MAP_INSTANTIATE(Real);

} // namespace libMesh
