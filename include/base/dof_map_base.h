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

#ifndef LIBMESH_DOF_MAP_BASE_H
#define LIBMESH_DOF_MAP_BASE_H

#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"
#include <vector>

namespace libMesh
{
class Variable;
class Elem;
class Node;

/**
 * This base class provides a minimal set of interfaces for satisfying user requests for
 * - variables that contribute degrees of freedom to this map
 * - distribution of degrees of freedom among ranks in parallel (e.g. MPI) simulations
 * - degrees of freedom associated with nodes and elements
 * This minimal set of interfaces is sufficient for providing information for use cases such
 * as building field splits for advanced preconditioning techniques
 *
 * Sub-classes of this base class may not be fully-featured stand-alone degree of freedom maps.
 * For instance a sub-class may not be constructible without information from an already constructed
 * fully-featured degree of freedom map. As an example, a degree of freedom map for a statically
 * condensed system of equations will require *a priori* and flow naturally from a degree of freedom
 * map for the full non-condensed system. In this example, the partitioning of the statically
 * condensed degrees of freedom may logically follow the partitioning in the global system, e.g. a
 * face degree of freedom in the global system that is assigned to process i by the full degree of
 * freedom map will retain its rank assignment in the statically condensed degree of freedom map.
 * Because the statically condensed dof map can lazily use the information from the full dof map,
 * it does not require its own partitioning capabilities. Consequently, partitioning interfaces
 * are not present in this base interface
 */
class DofMapBase : public ParallelObject
{
public:
  DofMapBase(const Parallel::Communicator & comm);

  /**
   * \returns The number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */
  virtual unsigned int n_variables() const = 0;

  /**
   * \returns The variable description object for variable \p c.
   */
  virtual const Variable & variable(const unsigned int c) const = 0;

  /**
   * \returns The first dof index that is local to partition \p proc.
   */
  dof_id_type first_dof(const processor_id_type proc) const;

  dof_id_type first_dof() const { return this->first_dof(this->processor_id()); }

  /**
   * \returns The first dof index that is after all indices local to
   * processor \p proc.
   *
   * Analogous to the end() member function of STL containers.
   */
  dof_id_type end_dof(const processor_id_type proc) const;

  dof_id_type end_dof() const { return this->end_dof(this->processor_id()); }

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element.  For one variable, and potentially for a
   * non-default element p refinement level
   */
  virtual void dof_indices(const Elem * const elem,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn,
                           int p_level = -12345) const = 0;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the \p node, for one variable \p vn.
   */
  virtual void dof_indices(const Node * const node,
                           std::vector<dof_id_type> & di,
                           const unsigned int vn) const = 0;
  /**
   * \returns The total number of degrees of freedom in the problem.
   */
  dof_id_type n_dofs() const { return _n_dfs; }

  /**
   * \returns The number of degrees of freedom on partition \p proc.
   */
  dof_id_type n_dofs_on_processor(const processor_id_type proc) const;

  /**
   * \returns The number of degrees of freedom on this processor.
   */
  dof_id_type n_local_dofs() const { return this->n_dofs_on_processor(this->processor_id()); }

  virtual void clear();

  /**
   * \returns The total number of degrees of freedom on old_dof_objects
   */
  dof_id_type n_old_dofs() const { return _n_old_dfs; }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * \returns The first old dof index that is local to partition \p proc.
   */
  dof_id_type first_old_dof(const processor_id_type proc) const;

  dof_id_type first_old_dof() const { return this->first_old_dof(this->processor_id()); }

  /**
   * \returns The first old dof index that is after all indices local
   * to processor \p proc.
   *
   * Analogous to the end() member function of STL containers.
   */
  dof_id_type end_old_dof(const processor_id_type proc) const;

  dof_id_type end_old_dof() const { return this->end_old_dof(this->processor_id()); }
#endif // LIBMESH_ENABLE_AMR

protected:
  /**
   * compute the key degree of freedom information given the local number of degrees of freedom on
   * this process
   * \returns The total number of DOFs for the System, summed across all procs.
   */
  std::size_t compute_dof_info(dof_id_type n_local_dofs);

  /**
   * First DOF index on processor \p p.
   */
  std::vector<dof_id_type> _first_df;

  /**
   * Last DOF index (plus 1) on processor \p p.
   */
  std::vector<dof_id_type> _end_df;

  /**
   * Total number of degrees of freedom.
   */
  dof_id_type _n_dfs;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Total number of degrees of freedom on old dof objects
   */
  dof_id_type _n_old_dfs;

  /**
   * First old DOF index on processor \p p.
   */
  std::vector<dof_id_type> _first_old_df;

  /**
   * Last old DOF index (plus 1) on processor \p p.
   */
  std::vector<dof_id_type> _end_old_df;
#endif
};

inline dof_id_type DofMapBase::first_dof(const processor_id_type proc) const
{
  libmesh_assert_less(proc, _first_df.size());
  return _first_df[proc];
}

inline dof_id_type DofMapBase::end_dof(const processor_id_type proc) const
{
  libmesh_assert_less(proc, _end_df.size());
  return _end_df[proc];
}

inline dof_id_type DofMapBase::n_dofs_on_processor(const processor_id_type proc) const
{
  libmesh_assert_less(proc, _first_df.size());
  return cast_int<dof_id_type>(_end_df[proc] - _first_df[proc]);
}

inline dof_id_type DofMapBase::first_old_dof(const processor_id_type proc) const
{
  libmesh_assert_less(proc, _first_old_df.size());
  return _first_old_df[proc];
}

inline dof_id_type DofMapBase::end_old_dof(const processor_id_type proc) const
{
  libmesh_assert_less(proc, _end_old_df.size());
  return _end_old_df[proc];
}

}
#endif // LIBMESH_DOF_MAP_BASE_H
