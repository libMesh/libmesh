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
#include "libmesh/elem.h"
#include "libmesh/topology_map.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"

// C++ Includes
#include <limits>
#include <utility>

namespace libMesh
{

void TopologyMap::init(MeshBase & mesh)
{
  // This function must be run on all processors at once
  // for non-serial meshes
  if (!mesh.is_serial())
    libmesh_parallel_only(mesh.comm());

  LOG_SCOPE("init()", "TopologyMap");

  // Clear the old map
  _map.clear();

  this->fill(mesh);
}



void TopologyMap::add_node(const Node & mid_node,
                           const std::vector<std::pair<dof_id_type, dof_id_type> > & bracketing_nodes)
{
  const dof_id_type mid_node_id = mid_node.id();

  libmesh_assert_not_equal_to(mid_node_id, DofObject::invalid_id);

  for (std::size_t i=0; i != bracketing_nodes.size(); ++i)
    {
      const dof_id_type id1 = bracketing_nodes[i].first;
      const dof_id_type id2 = bracketing_nodes[i].second;
      const dof_id_type lower_id = std::min(id1, id2);
      const dof_id_type upper_id = std::max(id1, id2);

      // We should never be inserting inconsistent data
#ifndef NDEBUG
      map_type::iterator it =
        _map.find(std::make_pair(lower_id, upper_id));

      if (it != _map.end())
        libmesh_assert_equal_to (it->second, mid_node_id);
#endif

      this->_map.insert(std::make_pair(std::make_pair(lower_id, upper_id),
                                       mid_node_id));

    }
}


dof_id_type TopologyMap::find(const std::vector<std::pair<dof_id_type, dof_id_type> > & bracketing_nodes) const
{
  dof_id_type new_node_id = DofObject::invalid_id;

  for (std::size_t i = 0; i != bracketing_nodes.size(); ++i)
    {
      const dof_id_type lower_id = std::min(bracketing_nodes[i].first,
                                            bracketing_nodes[i].second);
      const dof_id_type upper_id = std::max(bracketing_nodes[i].first,
                                            bracketing_nodes[i].second);

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



#ifdef LIBMESH_ENABLE_AMR

void TopologyMap::fill(const MeshBase & mesh)
{
  // Populate the nodes map
  MeshBase::const_element_iterator
    it = mesh.elements_begin(),
    end = mesh.elements_end();
  for (; it != end; ++it)
    {
      const Elem * elem = *it;

      // We only need to add nodes which might be added during mesh
      // refinement; this means they need to be child nodes.
      if (!elem->has_children())
        continue;

      for (unsigned int c = 0; c != elem->n_children(); ++c)
        {
          if (elem->child_ptr(c)->is_remote())
            continue;

          for (unsigned int n = 0; n != elem->n_nodes_in_child(c); ++n)
            {
              const std::vector<std::pair<dof_id_type, dof_id_type> >
                bracketing_nodes = elem->bracketing_nodes(c,n);

              this->add_node(elem->child_ptr(c)->node_ref(n),
                             bracketing_nodes);
            }
        }
    }
}

#else

// no-op without AMR
void TopologyMap::fill(const MeshBase &) {}

#endif



} // namespace libMesh
