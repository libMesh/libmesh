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

#ifndef LIBMESH_REDUCED_DOF_MAP_H
#define LIBMESH_REDUCED_DOF_MAP_H

#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map_base.h"
#include "libmesh/variable.h"
#include "libmesh/sparsity_pattern.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace libMesh
{
class MeshBase;
class System;
class DofMap;
class Elem;
class Node;

/**
 * A degree of freedom map for a reduced space: a subset of the degrees of freedom of a full \p
 * DofMap, renumbered contiguously. Each retained degree of freedom keeps the rank that owns it in
 * the full space, so the reduced numbering on each rank preserves the full numbering order.
 *
 * Subclasses select the retained degrees of freedom through \p retain_field_dof() and \p
 * retain_scalar_dof(). Static condensation, for instance, retains the degrees of freedom on
 * element boundaries, and constraint elimination retains the unconstrained ones. The constraining
 * degrees of freedom of every constrained degree of freedom on an element are retained as well,
 * whether or not the constrained one is, because constrained element matrices and vectors couple
 * them.
 *
 * This class provides the reduced numbering (\p reduced_dof()), the reduced degrees of freedom of
 * each local element and variable (\p dof_indices()), the variables present in the reduced space
 * along with a stand-in \p System for them, and the sparsity pattern of matrices in the reduced
 * space. This is the information needed to assemble reduced matrices and to build field splits
 * of the reduced space.
 */
class ReducedDofMap : public DofMapBase
{
public:
  virtual ~ReducedDofMap();

  /**
   * Build the reduced space from the current state of the full degree of freedom map
   */
  virtual void reinit();

  virtual void clear() override;

  /**
   * Whether \p reinit() has built the reduced space
   */
  bool initialized() const { return _initialized; }

  virtual unsigned int n_variables() const override;

  virtual const Variable & variable(const unsigned int c) const override;

  /**
   * Fills \p di with the reduced degree of freedom indices of variable \p vn on the active local
   * element \p elem, including the constraining degrees of freedom of its constrained degrees of
   * freedom
   */
  virtual void dof_indices(const Elem * const elem,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn,
                           int p_level = -12345) const override;

  /**
   * Reduced degree of freedom maps are queried by element only
   */
  virtual void dof_indices(const Node * const node,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn) const override;

  /**
   * \returns The reduced index of \p full_dof. \p full_dof must be retained and must be either
   * local or used by a local element
   */
  dof_id_type reduced_dof(dof_id_type full_dof) const;

  /**
   * \returns The local retained degrees of freedom in the full numbering, sorted. The reduced
   * index of entry \p i is <tt>first_dof() + i</tt>
   */
  const std::vector<dof_id_type> & local_retained_dofs() const { return _local_retained_dofs; }

  /**
   * Nests this reduced space in \p parent's, or in the full space if \p parent is null. Only
   * degrees of freedom that \p parent retains can be retained, and the parent space is the space
   * this one reduces: its vectors are the ones operators on this space are prolonged to and
   * restricted from. Takes effect at the next \p reinit(), which must follow the parent's
   */
  void set_parent(const ReducedDofMap * parent) { _parent = parent; }

  /**
   * \returns The reduced space this one is nested in, or null for the full space
   */
  const ReducedDofMap * parent() const { return _parent; }

  /**
   * \returns The index of \p full_dof in the parent space, or \p DofObject::invalid_id if the
   * parent does not retain it. \p full_dof must be local or used by an active local element
   */
  dof_id_type parent_dof(dof_id_type full_dof) const;

  /**
   * \returns The number of degrees of freedom in the parent space, globally and locally, and the
   * parent index range of this rank
   */
  dof_id_type parent_n_dofs() const;
  dof_id_type parent_n_local_dofs() const;
  dof_id_type parent_first_dof() const;
  dof_id_type parent_end_dof() const;

  /**
   * \returns The local retained degrees of freedom in the parent numbering, sorted. Entry \p i is
   * the parent index of reduced index <tt>first_dof() + i</tt>
   */
  const std::vector<dof_id_type> & local_retained_parent_dofs() const
  { return _local_retained_parent_dofs; }

  /**
   * \returns The nonlocal reduced degrees of freedom used by active local elements, sorted: the
   * ghost entries of a ghosted vector on this space
   */
  const std::vector<dof_id_type> & send_list() const { return _send_list; }

  /**
   * \returns The ghost entries of a ghosted vector on the parent space
   */
  const std::vector<dof_id_type> & parent_send_list() const;

  /**
   * Whether field degree of freedom \p field_dof is retained: the parent retains it and \p
   * retain_field_dof() accepts it
   */
  bool retains_field_dof(const Elem & elem,
                         unsigned int node_num,
                         unsigned int var_num,
                         dof_id_type field_dof) const;

  /**
   * Whether SCALAR degree of freedom \p scalar_dof is retained: the parent retains it and \p
   * retain_scalar_dof() accepts it
   */
  bool retains_scalar_dof(dof_id_type scalar_dof) const;

  /**
   * \returns The reduced degree of freedom indices on the active local element with id \p
   * elem_id, by reduced variable number. Concatenated, these are in the order in which \p
   * retained_dof_added() saw them
   */
  const std::vector<std::vector<dof_id_type>> & elem_dof_indices(dof_id_type elem_id) const;

  /**
   * \returns The stand-in system that holds the variables of the reduced space
   */
  const System & reduced_system() const;

  /**
   * \returns The mesh
   */
  const MeshBase & mesh() const { return _mesh; }

  /**
   * \returns The sparsity pattern of matrices in the reduced space. Its rows are indexed by the
   * local full degrees of freedom; use \p reduced_n_nz() and \p reduced_n_oz() for preallocation
   */
  const SparsityPattern::Build & reduced_sparsity() const;

  /**
   * \returns The number of on-processor nonzeros in each local row of the reduced space
   */
  const std::vector<dof_id_type> & reduced_n_nz() const { return _reduced_nnz; }

  /**
   * \returns The number of off-processor nonzeros in each local row of the reduced space
   */
  const std::vector<dof_id_type> & reduced_n_oz() const { return _reduced_noz; }

  /**
   * \returns The reduced index of \p full_dof, or \p DofObject::invalid_id if \p full_dof is
   * not retained. \p full_dof must be local or used by an active local element
   */
  dof_id_type find_reduced_dof(dof_id_type full_dof) const;

  /**
   * This space's rule for retaining field degree of freedom \p field_dof of variable \p var_num on
   * \p elem, among those the parent retains. \p node_num is the element-local node holding it, or
   * \p invalid_uint for an element degree of freedom
   */
  virtual bool retain_field_dof(const Elem & elem,
                                unsigned int node_num,
                                unsigned int var_num,
                                dof_id_type field_dof) const = 0;

  /**
   * This space's rule for retaining degree of freedom \p scalar_dof of a SCALAR variable, among
   * those the parent retains
   */
  virtual bool retain_scalar_dof(dof_id_type /*scalar_dof*/) const { return true; }

protected:
  /**
   * @param reduced_system_name The name of the stand-in system holding the reduced variables
   */
  ReducedDofMap(MeshBase & mesh,
                System & system,
                const DofMap & dof_map,
                std::string reduced_system_name);

  /**
   * Called by \p reinit() each time \p full_dof is first retained on \p elem, in the order of the
   * element's reduced indices
   */
  virtual void retained_dof_added(const Elem & elem, dof_id_type full_dof);

  /**
   * Called by \p reinit() for each field degree of freedom on \p elem that the parent retains and
   * \p retain_field_dof() rejects
   */
  virtual void unretained_dof_found(const Elem & elem, dof_id_type full_dof);

  MeshBase & _mesh;
  System & _system;
  const DofMap & _dof_map;

private:
  /**
   * Retain \p full_dof on \p elem, along with its constraining degrees of freedom, recursively
   * @param elem_seen The full degrees of freedom already retained on \p elem
   * @param elem_var_dofs The retained full degrees of freedom of the current variable on \p elem
   */
  void retain_dof_and_constraining_dofs(const Elem * elem,
                                        dof_id_type full_dof,
                                        std::unordered_set<dof_id_type> & elem_seen,
                                        std::vector<dof_id_type> & elem_var_dofs);

  /**
   * Retain on \p elem the constraining degrees of freedom of \p full_dof, if it is constrained
   */
  void retain_constraining_dofs(const Elem * elem,
                                dof_id_type full_dof,
                                std::unordered_set<dof_id_type> & elem_seen,
                                std::vector<dof_id_type> & elem_var_dofs);

  /**
   * Record that \p full_dof is retained, on this rank or on its owner
   */
  void retain_dof(dof_id_type full_dof);

  /// The name of the stand-in reduced system
  const std::string _reduced_system_name;

  /// The reduced space this one is nested in, or null for the full space
  const ReducedDofMap * _parent;

  /// The local retained degrees of freedom in the parent numbering, sorted
  std::vector<dof_id_type> _local_retained_parent_dofs;

  /// The nonlocal reduced degrees of freedom used by active local elements, sorted
  std::vector<dof_id_type> _send_list;

  /// Whether \p reinit() has built the reduced space
  bool _initialized;

  /// The local retained degrees of freedom in the full numbering, sorted
  std::vector<dof_id_type> _local_retained_dofs;

  /// The reduced index of each local full degree of freedom, offset by the first local full degree
  /// of freedom. Unretained degrees of freedom map to \p DofObject::invalid_id
  std::vector<dof_id_type> _local_full_to_reduced;

  /// The reduced indices of the nonlocal degrees of freedom on local elements, \p
  /// DofObject::invalid_id for unretained ones
  std::unordered_map<dof_id_type, dof_id_type> _nonlocal_full_to_reduced;

  /// Local retained degrees of freedom, nonlocal retained ones by owner, and nonlocal unretained
  /// ones on local elements. Only used during \p reinit()
  std::unordered_set<dof_id_type> _local_retained_set;
  std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> _nonlocal_retained;
  std::unordered_set<dof_id_type> _nonlocal_unretained;

  /// The reduced degrees of freedom of each active local element, by reduced variable number
  std::unordered_map<dof_id_type, std::vector<std::vector<dof_id_type>>> _elem_reduced_dofs;

  /// The variables in the reduced space
  std::vector<Variable> _reduced_vars;

  /// A stand-in system to give the reduced variables a system number
  System * _reduced_system;

  /// The reduced space sparsity pattern
  std::unique_ptr<SparsityPattern::Build> _reduced_sp;

  /// Number of on-processor nonzeros in each local row of the reduced space
  std::vector<dof_id_type> _reduced_nnz;

  /// Number of off-processor nonzeros in each local row of the reduced space
  std::vector<dof_id_type> _reduced_noz;
};

inline const System & ReducedDofMap::reduced_system() const
{
  libmesh_assert(_reduced_system);
  return *_reduced_system;
}

inline const SparsityPattern::Build & ReducedDofMap::reduced_sparsity() const
{
  libmesh_assert(_reduced_sp);
  return *_reduced_sp;
}

}

#endif // LIBMESH_REDUCED_DOF_MAP_H
