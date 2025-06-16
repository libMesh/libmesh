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

#ifndef LIBMESH_STATIC_CONDENSATION_DOF_MAP_H
#define LIBMESH_STATIC_CONDENSATION_DOF_MAP_H

#include "libmesh/libmesh_config.h"

#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map_base.h"
#include "libmesh/variable.h"
#include "libmesh/sparsity_pattern.h"
#include "libmesh/utility.h"

#include <unordered_map>
#include <memory>
#include <vector>

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
class StaticCondensationDofMap : public DofMapBase
{
public:
  StaticCondensationDofMap(MeshBase & mesh, System & system, const DofMap & dof_map);
  virtual ~StaticCondensationDofMap();

  //
  // Unique methods
  //

  /**
   * Build the element global to local index maps
   */
  void reinit();

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

  virtual unsigned int n_variables() const override;

  virtual const Variable & variable(const unsigned int c) const override;

  virtual void dof_indices(const Elem * const elem,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn,
                           int p_level = -12345) const override;

  virtual void dof_indices(const Node * const node,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn) const override;

  /*
   * @returns The dummyish reduced system
   */
  const System & reduced_system() const;

  /**
   * Whether we are initialized
   */
  bool initialized() const { return _sc_is_initialized; }

  virtual void clear() override;

  /**
   * Retrieve the dof index in the reduced system corresponding to the provided dof index in the
   * full system. This dof index should be a dof index involved in constraints
   */
  dof_id_type get_reduced_from_global_constraint_dof(dof_id_type full_dof) const;

private:
  /**
   * Add an uncondensed dof potentially along with constraining dofs which themselves must/will also
   * be uncondensed
   * @param full_dof_number The dof id in the full (uncondensed + condensed) system
   * @param involved_in_constraints Whether the \p full_dof_number is involved in constraints
   * @param uncondensed_local_to_global_map A map from uncondensed full dof numbering to local
   * (element) numbering for a given element
   * @param local_uncondensed_dofs_set Set of local uncondensed degrees of freedom (numbered in the
   * full system)
   * @param nonlocal_uncondensed_dofs Map from processor IDs to sets of nonlocal uncondensed degrees
   * of freedom (numbered in the full system)
   * @param elem_uncondensed_dofs Used for temporary storage of uncondensed degrees of freedom
   * active on an element
   * @param uncondensed_local_dof_number Counter for number of uncondensed dofs on an element
   * @param constraint_dofs Set of degrees of freedom numbered in the full system involved in
   * constraints
   */
  void add_uncondensed_dof_plus_constraint_dofs(
      dof_id_type full_dof_number,
      bool involved_in_constraints,
      std::unordered_map<dof_id_type, dof_id_type> & uncondensed_global_to_local_map,
      std::unordered_set<dof_id_type> & local_uncondensed_dofs_set,
      std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> &
          nonlocal_uncondensed_dofs,
      std::vector<dof_id_type> & elem_uncondensed_dofs,
      dof_id_type & uncondensed_local_dof_number,
      std::unordered_set<dof_id_type> & constraint_dofs);

  /**
   * Add an uncondensed dof
   * @param full_dof_number The dof id in the full (uncondensed + condensed) system
   * @param involved_in_constraints Whether the \p full_dof_number is involved in constraints
   * @param uncondensed_local_to_global_map A map from uncondensed full dof numbering to local
   * (element) numbering for a given element
   * @param local_uncondensed_dofs_set Set of local uncondensed degrees of freedom (numbered in the
   * full system)
   * @param nonlocal_uncondensed_dofs Map from processor IDs to sets of nonlocal uncondensed degrees
   * of freedom (numbered in the full system)
   * @param elem_uncondensed_dofs Used for temporary storage of uncondensed degrees of freedom
   * active on an element
   * @param uncondensed_local_dof_number Counter for number of uncondensed dofs on an element
   * @param constraint_dofs Set of degrees of freedom numbered in the full system involved in
   * constraints
   */
  void add_uncondensed_dof(
      dof_id_type full_dof_number,
      bool involved_in_constraints,
      std::unordered_map<dof_id_type, dof_id_type> & uncondensed_global_to_local_map,
      std::unordered_set<dof_id_type> & local_uncondensed_dofs_set,
      std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> &
          nonlocal_uncondensed_dofs,
      std::vector<dof_id_type> & elem_uncondensed_dofs,
      dof_id_type & uncondensed_local_dof_number,
      std::unordered_set<dof_id_type> & constraint_dofs);

  /**
   * Data stored on a per-element basis used to compute element Schur complements and their
   * applications to vectors
   */
  struct DofData
  {
    /// The uncondensed degrees of freedom with global numbering corresponding to the the \emph reduced
    /// system. Note that initially this will actually hold the indices corresponding to the fully
    /// sized problem, but we will swap it out by the time we are done initializing
    std::vector<std::vector<dof_id_type>> reduced_space_indices;

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
  };

  /// Friend the static condensation matrix class so they can inspect our degree of freedom data
  friend class StaticCondensation;

  /// A map from element ID to Schur complement data
  std::unordered_map<dof_id_type, DofData> _elem_to_dof_data;

  /// All the uncondensed degrees of freedom (numbered in the "full" uncondensed + condensed
  /// space). This data member is used for creating subvectors corresponding to only uncondensed
  /// dofs
  std::vector<dof_id_type> _local_uncondensed_dofs;

  MeshBase & _mesh;
  System & _system;
  const DofMap & _dof_map;

  /// Variables for which we will keep all dofs
  std::unordered_set<unsigned int> _uncondensed_vars;

  /// Whether our object has been initialized
  bool _sc_is_initialized;

  /// The variables in the reduced system
  std::vector<Variable> _reduced_vars;

  /// A dummyish system to help with DofObjects
  System * _reduced_system;

  /// Owned storage of the reduced system sparsity pattern
  std::unique_ptr<SparsityPattern::Build> _reduced_sp;

  /// Number of on-diagonal nonzeros per row in the reduced system
  std::vector<dof_id_type> _reduced_nnz;

  /// Number of off-diagonal nonzeros per row in the reduced system
  std::vector<dof_id_type> _reduced_noz;

  /// A small map from full system degrees of freedom to reduced/condensed system degrees of freedom
  /// involved in constraints. This may be leveraged by \p StaticCondensation::set during \p
  /// DofMap::enforce_constraints_on_jacobian at the very end of Jacobian evaluation
  std::unordered_map<dof_id_type, dof_id_type> _full_to_reduced_constraint_dofs;
};

inline void
StaticCondensationDofMap::dont_condense_vars(const std::unordered_set<unsigned int> & vars)
{
  _uncondensed_vars.insert(vars.begin(), vars.end());
}

inline const System & StaticCondensationDofMap::reduced_system() const
{
  libmesh_assert(_reduced_system);
  return *_reduced_system;
}

inline dof_id_type
StaticCondensationDofMap::get_reduced_from_global_constraint_dof(const dof_id_type full_dof) const
{
  return libmesh_map_find(_full_to_reduced_constraint_dofs, full_dof);
}

}

#endif // LIBMESH_STATIC_CONDENSATION_DOF_MAP_H
