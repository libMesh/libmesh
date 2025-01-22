// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// libMesh includes
#include "libmesh/sides_to_elem_map.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h" // LOG_SCOPE
#include "libmesh/mesh_base.h"
#include "libmesh/utility.h" // libmesh_map_find()

namespace libMesh
{

namespace MeshTools
{

SidesToElemMap::SidesToElemMap() = default;
SidesToElemMap::~SidesToElemMap() = default;

/**
 * Static build function
 */
SidesToElemMap SidesToElemMap::build(const MeshBase & mesh)
{
  LOG_SCOPE("build()", "SidesToElemMap");

  // Eventual return value
  SidesToElemMap ret;

  // Temporary storage of sorted list of vertex ids on side "s" of
  // "elem".  We use only vertices since that is sufficient to
  // uniquely identify an Elem side. We only create this once and then
  // reuse it for each Side in the mesh to avoid repeated small memory
  // allocations.
  std::vector<dof_id_type> sorted_vertex_ids;

  // Flag which is used to check for hanging Nodes, a case not
  // currently handled by the SidesToElemMap object.
  unsigned int active_level_seen = libMesh::invalid_uint;

  for (const auto & elem : mesh.element_ptr_range())
  {
    // We currently don't handle meshes with hanging nodes. That is,
    // if there is a refined element with a coarser neighbor, the
    // SidesToElemMap won't contain any neighbor information (in
    // either direction) about this. For now, we throw an error
    // in this scenario.
    if (elem->active())
      {
        if (active_level_seen == libMesh::invalid_uint)
          active_level_seen = elem->level();
        else
          libmesh_error_msg_if(
            active_level_seen != elem->level(),
            "SidesToElemMap does not currently support hanging nodes.");
      }

    for (auto s : elem->side_index_range())
    {
      ret.get_sorted_vertex_ids(elem, s, sorted_vertex_ids);

      // Get reference to (or create) the vector of Elem pointers
      // associated with this list of vertex ids, and store "elem" in it.
      ret._sides_to_elem_map[sorted_vertex_ids].push_back(elem);
    }
  }

  return ret;
}

std::pair<SidesToElemMap::ElemIter, SidesToElemMap::ElemIter>
SidesToElemMap::get_connected_elems(const Elem * elem, unsigned int side) const
{
  std::vector<dof_id_type> sorted_vertex_ids;
  this->get_sorted_vertex_ids(elem, side, sorted_vertex_ids);

  const auto & elem_vec = libmesh_map_find(_sides_to_elem_map, sorted_vertex_ids);
  return std::make_pair(elem_vec.begin(), elem_vec.end());
}

void
SidesToElemMap::get_sorted_vertex_ids(
  const Elem * elem,
  unsigned int side,
  std::vector<dof_id_type> & sorted_vertex_ids) const
{
  // Clear any prior data
  sorted_vertex_ids.clear();

  // Get all vertices for this side. Note: we can stop once we reach
  // the first non-vertex node, since all vertices are numbered first.
  for (const auto & n : elem->nodes_on_side(side))
  {
    if (elem->is_vertex(n))
      sorted_vertex_ids.push_back(elem->node_id(n));
    else
      break;
  }

  // Sort the vertex ids. The hash functions depend on the ordering
  // of the ids, so we want to be sure they are the same for every
  // Elem that touches this Side.
  std::sort(sorted_vertex_ids.begin(), sorted_vertex_ids.end());
}

} // namespace MeshTools

} // namespace libMesh

