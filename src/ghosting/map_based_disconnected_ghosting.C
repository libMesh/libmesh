#include "libmesh/map_based_disconnected_ghosting.h"

#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/simple_range.h"

#include <memory>

namespace libMesh
{

MapBasedDisconnectedGhosting::
MapBasedDisconnectedGhosting(MeshBase & mesh,
                               const DisconnectedMap & map) :
  GhostingFunctor(mesh),
  _map(map)
{
}

std::unique_ptr<GhostingFunctor>
MapBasedDisconnectedGhosting::clone () const
{
  // Create a new unique_ptr by copying this object
  return std::make_unique<MapBasedDisconnectedGhosting>(*this);
}

void
MapBasedDisconnectedGhosting::operator() (
  const MeshBase::const_element_iterator & range_begin,
  const MeshBase::const_element_iterator & range_end,
  processor_id_type p,
  map_type & coupled_elements)
{
  CouplingMatrix * nullcm = nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      // Add the element from the input range (if not local)
      if (elem->processor_id() != p)
        coupled_elements.emplace(elem, nullcm);

      // Look up this element's neighbor in the map
      auto it = _map.find(elem->id());
      if (it == _map.end())
        continue;

      // Get the neighbor's ID
      const dof_id_type neighbor_id = it->second;

      const Elem * neighbor_ptr = _mesh->query_elem_ptr(neighbor_id);

      // If the neighbor pointer is valid (meaning it's local or ghosted)
      // AND it's not on processor 'p', add it to the couplings.
      if (neighbor_ptr && neighbor_ptr->processor_id() != p)
        {
          coupled_elements.emplace(neighbor_ptr, nullcm);
        }
    }
}

} // namespace libMesh
