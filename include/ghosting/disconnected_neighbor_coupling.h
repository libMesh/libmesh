#ifndef LIBMESH_DISCONNECTED_NEIGHBOR_COUPLING_H
#define LIBMESH_DISCONNECTED_NEIGHBOR_COUPLING_H

#include "libmesh/ghosting_functor.h"
#include "libmesh/boundary_info.h"

#include <memory>

namespace libMesh
{

class Elem;
class PointLocatorBase;

/**
 * Ghosts elements that lie on
 * user-defined disconnected interfaces,
 * such as cohesive zone boundaries.
 */
class DisconnectedNeighborCoupling : public GhostingFunctor
{
public:
  using DisconnectedMap = std::unordered_map<const Elem *, const Elem *>;

  DisconnectedNeighborCoupling(const MeshBase & mesh) :
    GhostingFunctor(mesh),
    _dof_coupling(nullptr)
  {}

  virtual std::unique_ptr<GhostingFunctor> clone () const override
  { return std::make_unique<DisconnectedNeighborCoupling>(*this); }

  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

private:
  const CouplingMatrix * _dof_coupling;
};

} // namespace libMesh

#endif // LIBMESH_DISCONNECTED_NEIGHBOR_COUPLING_H
