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



#ifndef LIBMESH_TOPOLOGY_MAP_H
#define LIBMESH_TOPOLOGY_MAP_H

// Local Includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"

// C++ Includes
#include <unordered_map>
#include <functional> // std::hash
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;
template <typename> class NodeTempl;
typedef NodeTempl<Real> Node;

// Fix for STL laziness
struct myhash {
public:
  template <typename T1, typename T2>
  std::size_t operator()(const std::pair<T1, T2> & x) const
  {
    // recommendation from
    // http://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
    return 3 * std::hash<T1>()(x.first) + std::hash<T2>()(x.second);
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
 *
 * \author Roy Stogner
 * \date 2015
 * \brief Enables topology-based lookups of nodes.
 */
class TopologyMap
{
  // We need to supply our own hash function.
  typedef std::unordered_map<std::pair<dof_id_type, dof_id_type>, dof_id_type, myhash> map_type;
public:
  template <typename RealType>
  void init(MeshBaseTempl<RealType> &);

  void clear() { _map.clear(); }

  /**
   * Add a node to the map, between each pair of specified bracketing
   * nodes.
   */
  template <typename RealType>
  void add_node(const NodeTempl<RealType> & mid_node,
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
  template <typename RealType>
  void fill(const MeshBaseTempl<RealType> &);

private:
  map_type          _map;
};

} // namespace libMesh


#endif // LIBMESH_TOPOLOGY_MAP_H
