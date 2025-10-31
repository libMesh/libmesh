#include "libmesh/disconnected_neighbor_coupling.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/remote_elem.h"
#include "libmesh/periodic_boundaries.h"

namespace libMesh
{

void DisconnectedNeighborCoupling::operator()(
    const MeshBase::const_element_iterator & range_begin,
    const MeshBase::const_element_iterator & range_end,
    processor_id_type p,
    map_type & coupled_elements)
{
  libmesh_assert(_mesh);
  auto * db = _mesh->get_disconnected_boundaries();
  if (!db)
    return;

  std::unique_ptr<PointLocatorBase> point_locator = _mesh->sub_point_locator();

  // Primary ghosting: elements in range_begin...range_end
  for (const auto & elem : as_range(range_begin, range_end))
    if (elem->processor_id() != p)
      coupled_elements.emplace(elem, _dof_coupling);

  // Also ghost their disconnected neighbors
  for (const auto & elem : as_range(range_begin, range_end))
    for (auto s : elem->side_index_range())
    {
      for (const auto & [id, boundary_ptr] : *db)
      {
        if (!_mesh->get_boundary_info().has_boundary_id(elem, s, id))
          continue;

        unsigned int neigh_side = invalid_uint;
        const Elem * neigh =
            db->neighbor(id, *point_locator, elem, s, &neigh_side);

        if (!neigh || neigh == remote_elem)
          continue;

        if (neigh->processor_id() != p)
          coupled_elements.emplace(neigh, _dof_coupling);
      }
    }
}


} // namespace libMesh
