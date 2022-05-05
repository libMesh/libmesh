// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TOPOLOGY_MAP_H
#define LIBMESH_TOPOLOGY_MAP_H

// Local Includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/hashing.h"

// C++ Includes
#include <unordered_map>
#include <vector>

namespace libMesh
{

// Forward Declarations
class Elem;
class MeshBase;
class Node;

/**
 * Data structures that enable topology-based lookups of nodes created
 * by mesh refinement.
 *
 * The key is a pair of node ids for two nodes bracketing the new
 * node, sorted lowest id first.
 *
 * A node created in the middle of a cell's quad face will be the
 * value of two keys, one for each node pair bracketing it.
 *
 * For efficiency we will use a hashed map if it is available,
 * otherwise a regular map.
 *
 * \author Roy Stogner
 * \date 2015
 * \brief Enables topology-based lookups of nodes.
 */
class TopologyMap
{
  // We need to supply our own hash function.
  typedef std::unordered_map<std::pair<dof_id_type, dof_id_type>, dof_id_type, libMesh::hash> map_type;
public:
  void init(MeshBase &);

  void clear() { _map.clear(); }

  /**
   * Add a node to the map, between each pair of specified bracketing
   * nodes.
   */
  void add_node(const Node & mid_node,
                const std::vector<
                std::pair<dof_id_type, dof_id_type>> &
                bracketing_nodes);

  bool empty() const { return _map.empty(); }

  dof_id_type find(dof_id_type bracket_node1,
                   dof_id_type bracket_node2) const;

  dof_id_type find(const std::vector<
                   std::pair<dof_id_type, dof_id_type>> &
                   bracketing_nodes) const;

protected:
  void fill(const MeshBase &);

private:
  map_type          _map;
};

} // namespace libMesh


#endif // LIBMESH_TOPOLOGY_MAP_H
