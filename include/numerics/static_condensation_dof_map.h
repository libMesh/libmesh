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

#ifndef LIBMESH_STATIC_CONDENSATION_DOF_MAP_H
#define LIBMESH_STATIC_CONDENSATION_DOF_MAP_H

#include "libmesh/libmesh_config.h"

#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/reduced_dof_map.h"

#include <unordered_map>
#include <unordered_set>

namespace libMesh
{
class MeshBase;
class System;
class DofMap;
class Elem;

/**
 * A class holding degree of freedom information pertinent to static condensation.
 * Static condensation is a process for reducing the number of unknowns in a
 * linear(ized) system. It is often used in conjunction with a hybridized discontinuous Galerkin
 * (HDG) method to reduce the number of unknowns in the linear system to something more on the order
 * of a continuous Galerkin (CG) method as opposed to a discontinuous Galerkin method. In the HDG
 * example, degrees of freedom on the element interiors are condensed out and only the degrees of
 * freedom on element faces are retained in the condensed system. Static condensation is not only
 * applicable to HDG methods, however. One could also use static condensation with a high order CG
 * method, again removing internal degrees of freedom and leaving those on element boundaries.
 *
 * Users may query this class for condensed space degree of freedom information associated with
 * elements. Queries with nodes are not supported. This class also holds information for mapping
 * degree of freedom indices in the global, uncondensed space to local element degree of freedom
 * numbers. This information is essential in order to build the element Schur complement matrices
 * that are used to map global system vector and matrix data to the condensed system vector and
 * matrix data
 */
class StaticCondensationDofMap : public ReducedDofMap
{
public:
  StaticCondensationDofMap(MeshBase & mesh, System & system, const DofMap & dof_map);
  virtual ~StaticCondensationDofMap();

  /**
   * Build the reduced space and the element global to local index maps
   */
  virtual void reinit() override;

  /**
   * Add \p vars to the list of variables not to condense. This can be useful when some variable's
   * equation is discretized with a DG method or if including the variable in the condensed block
   * diagonal would result in it being singular
   */
  void dont_condense_vars(const std::unordered_set<unsigned int> & vars);

  /**
   * @returns our list of variables for whom we do not condense out any dofs
   */
  const std::unordered_set<unsigned int> & uncondensed_vars() const { return _uncondensed_vars; }

  virtual void clear() override;

  /**
   * Degrees of freedom on element boundaries and degrees of freedom of uncondensed variables are
   * retained
   */
  virtual bool retain_field_dof(const Elem & elem,
                                unsigned int node_num,
                                unsigned int var_num,
                                dof_id_type field_dof) const override;

protected:
  virtual void retained_dof_added(const Elem & elem, dof_id_type full_dof) override;

  virtual void unretained_dof_found(const Elem & elem, dof_id_type full_dof) override;

private:
  /**
   * Data stored on a per-element basis used to compute element Schur complements and their
   * applications to vectors
   */
  struct DofData
  {
    /// A map from the global degree of freedom number for the full system (condensed + uncondensed)
    /// to an element local number. If this map is queried with a condensed dof, nothing will be
    /// found. The size of this container will be the number of uncondensed degrees of freedom whose
    /// basis functions are nonzero on the element
    std::unordered_map<dof_id_type, dof_id_type> uncondensed_global_to_local_map;
    /// A map from the global degree of freedom number for the full system (condensed + uncondensed)
    /// to an element local number. If this map is queried with an uncondensed dof, nothing will be
    /// found. The size of this container will be the number of condensed degrees of freedom whose
    /// basis functions are nonzero on the element
    std::unordered_map<dof_id_type, dof_id_type> condensed_global_to_local_map;
    /// The parent-space indices of the uncondensed and condensed degrees of freedom, by element
    /// local number
    std::vector<dof_id_type> uncondensed_parent_dofs;
    std::vector<dof_id_type> condensed_parent_dofs;
  };

  /// Friend the static condensation matrix class so they can inspect our degree of freedom data
  friend class StaticCondensation;

  /// A map from element ID to Schur complement data
  std::unordered_map<dof_id_type, DofData> _elem_to_dof_data;

  /// Variables for which we will keep all dofs
  std::unordered_set<unsigned int> _uncondensed_vars;
};

inline void
StaticCondensationDofMap::dont_condense_vars(const std::unordered_set<unsigned int> & vars)
{
  _uncondensed_vars.insert(vars.begin(), vars.end());
}
}

#endif // LIBMESH_STATIC_CONDENSATION_DOF_MAP_H
