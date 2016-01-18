// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"

// C++ Includes   -----------------------------------
#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_HASH
#include <vector>

namespace libMesh
{

// Forward Declarations -----------------------------
class Elem;
class MeshBase;
class Node;

// Fix for STL laziness
struct myhash {
public:
  template <typename T1, typename T2>
  std::size_t operator()(const std::pair<T1, T2> & x) const
  {
    // recommendation from
    // http://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
    return 3 * LIBMESH_BEST_HASH<T1>()(x.first) + LIBMESH_BEST_HASH<T2>()(x.second);
  }
};

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
 */
class TopologyMap
{
  // We need to supply our own hash function if we're hashing
#if defined(LIBMESH_HAVE_STD_UNORDERED_MAP) ||  \
  defined(LIBMESH_HAVE_TR1_UNORDERED_MAP) ||    \
  defined(LIBMESH_HAVE_EXT_HASH_MAP) ||         \
  defined(LIBMESH_HAVE_HASH_MAP)
#  define MYHASH ,myhash
#else
#  define MYHASH
#endif

  typedef LIBMESH_BEST_UNORDERED_MAP<std::pair<dof_id_type, dof_id_type>,
                                     dof_id_type MYHASH> map_type;
public:
  void init(MeshBase &);

  void clear() { _map.clear(); }

  /**
   * Add a node to the map, between each pair of specified bracketing
   * nodes.
   */
  void add_node(const Node & mid_node,
                const std::vector<
                std::pair<dof_id_type, dof_id_type> > &
                bracketing_nodes);

  bool empty() const { return _map.empty(); }

  dof_id_type find(dof_id_type bracket_node1,
                   dof_id_type bracket_node2) const;

  dof_id_type find(const std::vector<
                   std::pair<dof_id_type, dof_id_type> > &
                   bracketing_nodes) const;

protected:
  void fill(const MeshBase &);

private:
  map_type          _map;
};

} // namespace libMesh


#endif // LIBMESH_TOPOLOGY_MAP_H
