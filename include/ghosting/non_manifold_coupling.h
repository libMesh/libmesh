// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_NON_MANIFOLD_COUPLING_H
#define LIBMESH_NON_MANIFOLD_COUPLING_H

// libMesh Includes
#include "libmesh/ghosting_functor.h" // base class
#include "libmesh/sides_to_elem_map.h" // MeshTools::SidesToElemMap

// C++ includes
#include <memory>

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * A GhostingFunctor subclass that uses a SidesToElemMap object to
 * determine when Elem DOFs should be ghosted in non-manifold meshes.
 * A basic use case for this functor is "truss" meshes where multiple
 * 1D elements meet at a single node, but there are also 2D
 * generalizations of practical importance, such as "shell stiffener"
 * models where three or more quadrilateral elements meet at the same
 * side. This functor can be used to avoid DOF ghosting issues related
 * to non-manifold meshes which are due to the assumption in libMesh
 * that there is at most one neighbor per Elem side.
 * Notes on the GhostingFunctor-derived classes:
 *
 * * For the GhostingFunctor::operator() overrides, the input range
 *   is saying, "we already need these Elems' DOFs, so what other
 *   Elems' DOFs do we need?"  This is why we always see the following
 *   logic in the GhostingFunctor::operator() overrides:
 *
 *      if (elem->processor_id() != p)
 *        coupled_elements.emplace(elem, _dof_coupling);
 *
 * * The idea is that the "output" range that we are creating should always
 *   include the entire "non-local" part of the input range. We don't need to
 *   include members of the input range that are already owned by processor p,
 *   see e.g. the API docs which say:
 *   "Don't bother to return any results which already have processor_id p."
 *
 * \author John W. Peterson
 * \date 2024
 */
class NonManifoldGhostingFunctor : public GhostingFunctor
{
public:
  /**
   * Constructor. Calls base class constructor.
   */
  NonManifoldGhostingFunctor(const MeshBase & mesh);

  /**
   * Special functions.
   * These are all defaulted for now until there is some reason to
   * delete/customize one of them.
   */
  NonManifoldGhostingFunctor (const NonManifoldGhostingFunctor &) = default;
  NonManifoldGhostingFunctor (NonManifoldGhostingFunctor &&) = default;
  NonManifoldGhostingFunctor & operator= (const NonManifoldGhostingFunctor &) = default;
  NonManifoldGhostingFunctor & operator= (NonManifoldGhostingFunctor &&) = default;
  virtual ~NonManifoldGhostingFunctor() = default;

  /**
   * clone() just calls the copy ctor.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const override;

  /**
   * If the Mesh changes, we'll need to clear the SidesToElemMap and
   * recreate it, since all the neighbor information (standard and
   * non-manifold) might have changed.
   */
  virtual void mesh_reinit () override;

  /**
   * When the Mesh is redistributed, we may need to update the
   * information in the SidesToElemMap by removing pointers to
   * non-local Elems. Currently calls mesh_reinit().
   */
  virtual void redistribute () override;

  /**
   * When remote elements are deleted, we should remove them from the
   * information in the SidesToElemMap. Currently calls mesh_reinit().
   */
  virtual void delete_remote_elements() override;

  /**
   * For the specified range of active elements, query the SidesToElemMap
   * and add _all_ side neighbors to the list of coupled_elements, ignoring
   * those on processor p.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;
private:

  /**
   * Quickly-searchable map from (elem, side) pair to list of connected Elems.
   */
  MeshTools::SidesToElemMap _stem;
};

} // namespace libMesh

#endif // LIBMESH_NON_MANIFOLD_COUPLING_H
