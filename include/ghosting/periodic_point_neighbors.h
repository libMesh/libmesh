// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PERIODIC_POINT_NEIGHBORS_H
#define LIBMESH_PERIODIC_POINT_NEIGHBORS_H

// Local Includes
#include "libmesh/ghosting_functor.h"

// C++ Includes
#include <memory>

namespace libMesh
{
#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;

/**
 * This class ghosts the point neighbors across periodic boundaries: the active elements on a
 * periodic boundary that touch the periodic image of a local element. Periodic constraint
 * computation compares every element at a periodic vertex and at its image, so it needs these
 * elements.
 *
 * Each DofMap adds one of these to its mesh for its own periodic boundaries.
 */
class PeriodicPointNeighbors : public GhostingFunctor
{
public:
  PeriodicPointNeighbors(const MeshBase & mesh) :
      GhostingFunctor(mesh),
      _periodic_bcs(nullptr)
    {}

  /**
   * The copy ghosts nothing until it receives periodic boundaries. The other's boundaries belong
   * to a DofMap that may be unrelated to any DofMap that depends on this object's ghosting.
   */
  PeriodicPointNeighbors(const PeriodicPointNeighbors & other) :
      GhostingFunctor(other),
      _periodic_bcs(nullptr)
    {}

  /**
   * A clone() is needed because GhostingFunctor can not be shared between
   * different meshes. The operations in  GhostingFunctor are mesh dependent.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const override
  { return std::make_unique<PeriodicPointNeighbors>(*this); }

  void set_periodic_boundaries(const PeriodicBoundaries * periodic_bcs) override
  { _periodic_bcs = periodic_bcs; }

  /**
   * If we have periodic boundaries, then we'll need the mesh to have
   * an updated point locator whenever we're about to query them.
   */
  virtual void mesh_reinit () override;

  virtual void redistribute () override
  { this->mesh_reinit(); }

  virtual void delete_remote_elements() override
  { this->mesh_reinit(); }

  /**
   * For the specified range of active elements, find their point neighbors across periodic
   * boundaries, ignoring those on processor p.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

private:
  const PeriodicBoundaries * _periodic_bcs;
};
#endif // LIBMESH_ENABLE_PERIODIC

} // namespace libMesh

#endif // LIBMESH_PERIODIC_POINT_NEIGHBORS_H
