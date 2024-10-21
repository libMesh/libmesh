// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

  for (const auto & elem : mesh.element_ptr_range())
    for (auto s : elem->side_index_range())
    {
      ret.get_sorted_vertex_ids(elem, s, sorted_vertex_ids);

      // Get reference to (or create) the vector of Elem pointers
      // associated with this list of vertex ids, and store "elem" in it.
      ret._sides_to_elem_map[sorted_vertex_ids].push_back(elem);
    }

  return ret;
}

const std::vector<const Elem *> &
SidesToElemMap::get_connected_elems(const Elem * elem, unsigned int side) const
{
  std::vector<dof_id_type> sorted_vertex_ids;
  this->get_sorted_vertex_ids(elem, side, sorted_vertex_ids);

  return libmesh_map_find(_sides_to_elem_map, sorted_vertex_ids);
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


void SidesToElemMap::print_statistics() const
{
  // An N-junction is defined as a Side which is connected to N elements.
  // Here we count number of occurrences of 2-, 3-, ... N-junctions.
  std::map<unsigned int, unsigned int> junction_counts;
  for (const auto & pr : _sides_to_elem_map)
    junction_counts[pr.second.size()]++;

  // Debugging:
  // Verify that we have the following number of N-junctions
  // (sides with N attached Elems)
  // N = 1, 964
  // N = 2, 7030
  // N = 3, 716
  // N = 4, 232
  // N = 5, 4
  // Notes:
  // * There are no N=0 sides since that would be a Side not attached to any Elems
  // * N=1 means it is a boundary side
  // * N=2 is the most common case where two Elems meet in a side
  // libMesh::out << "Number of sides with N attached Elems:" << std::endl;
  for (const auto & [jct_size, cnt] : junction_counts)
    libMesh::out << "N = " << jct_size << ", " << cnt << std::endl;
}

void SidesToElemMap::print_info() const
{
  libMesh::out << "Sets of connected Elems:" << std::endl;
  for (const auto & pr : _sides_to_elem_map)
  {
    for (const auto & elem : pr.second)
      libMesh::out << elem->id() << " ";
    libMesh::out << std::endl;
  }
}

std::vector<std::pair<dof_id_type, dof_id_type>>
SidesToElemMap::search_for_separated_elems() const
{
  std::vector<std::pair<dof_id_type, dof_id_type>> ret;

  // Search for pairs of Elems that:
  // 1.) Share a Side shared by > 1 other neighbors
  // 2.) Are not "standard" face neighbors
  // 3.) Are partitioned onto different procs
  for (const auto & pr : _sides_to_elem_map)
  {
    const auto & elem_vec = pr.second;

    // Only care about Sides with 3 or more neighbors
    if (elem_vec.size() < 3)
      continue;

    // Convenient pointer to first Elem in the list
    const Elem * elem0 = elem_vec[0];
    processor_id_type elem0_proc = elem0->processor_id();

    // An in-place Lambda that checks if elem_vec[0]'s processor_id is
    // different from any other Elem in elem_vec, and which aren't
    // already "standard" face neighbors, since those are already
    // ghosted.  Note that it doesn't try to find every pair that
    // differs, just a candidate pair.
    auto e = [&]()
      {
        for (unsigned int e=1; e<elem_vec.size(); ++e)
          if (elem_vec[e]->processor_id() != elem0_proc &&
              elem_vec[e]->which_neighbor_am_i(elem0) == libMesh::invalid_uint)
            return e;

        // Didn't find a pair owned by different processors that
        // weren't standard face neighbors
        return libMesh::invalid_uint;
      }();

    if (e != libMesh::invalid_uint)
    {
      // libMesh::out << "Elems " << elem0->id() << " and " << elem_vec[e]->id()
      //              << " share a Side but are on different processors and are not standard face neighbors." << std::endl;
      ret.emplace_back(elem0->id(), elem_vec[e]->id());
    }
  }

  return ret;
}

} // namespace MeshTools

} // namespace libMesh

