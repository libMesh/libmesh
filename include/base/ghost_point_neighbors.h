// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_GHOST_POINT_NEIGHBORS_H
#define LIBMESH_GHOST_POINT_NEIGHBORS_H

// Local Includes
#include "libmesh/ghosting_functor.h"
#include "libmesh/auto_ptr.h"

namespace libMesh
{
#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
#endif

/**
 * This class implements the original default geometry ghosting
 * requirements in libMesh: point neighbors on the same manifold, and
 * interior_parent elements.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
class GhostPointNeighbors : public GhostingFunctor
{
public:

  /**
   * Constructor.
   */
  GhostPointNeighbors(const MeshBase & mesh) :
      GhostingFunctor(mesh)
#ifdef LIBMESH_ENABLE_PERIODIC
      ,
      _periodic_bcs(nullptr)
#endif
    {}

  /**
   * Constructor.
   */
  GhostPointNeighbors(const GhostPointNeighbors & other) :
      GhostingFunctor(other)
#ifdef LIBMESH_ENABLE_PERIODIC
        ,
      // We do not simply want to copy over the other's periodic bcs because
      // they may very well correspond to periodic bcs from a \p DofMap entirely
      // unrelated to any \p DofMaps that depend on this objects ghosting
      _periodic_bcs(nullptr)
#endif
    {}

  /**
   * A clone() is needed because GhostingFunctor can not be shared between
   * different meshes. The operations in  GhostingFunctor are mesh dependent.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const override
  { return libmesh_make_unique<GhostPointNeighbors>(*this); }

#ifdef LIBMESH_ENABLE_PERIODIC
  // Set PeriodicBoundaries to couple
  void set_periodic_boundaries(const PeriodicBoundaries * periodic_bcs) override
  { _periodic_bcs = periodic_bcs; }
#endif

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
   * For the specified range of active elements, find their point
   * neighbors and interior_parent elements, ignoring those on
   * processor p.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

private:
#ifdef LIBMESH_ENABLE_PERIODIC
  const PeriodicBoundaries * _periodic_bcs;
#endif
};

} // namespace libMesh

#endif // LIBMESH_GHOST_POINT_NEIGHBORS_H
