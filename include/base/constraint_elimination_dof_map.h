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

#ifndef LIBMESH_CONSTRAINT_ELIMINATION_DOF_MAP_H
#define LIBMESH_CONSTRAINT_ELIMINATION_DOF_MAP_H

#include "libmesh/reduced_dof_map.h"

namespace libMesh
{

/**
 * The reduced space of the unconstrained degrees of freedom of a \p DofMap.
 *
 * Every \p DofMap constraint is affine: the full solution is \f$ x = C y + h \f$, where \f$ y \f$
 * holds the unconstrained degrees of freedom, \f$ C \f$ is the identity on them and holds the
 * constraint coefficients on constrained rows, and \f$ h \f$ holds the constraint values. The
 * reduced residual and Jacobian are \f$ C^T F(C y + h) \f$ and \f$ C^T J C \f$. Element matrices
 * and vectors constrained with \p DofMap::constrain_element_matrix() and its relatives are \f$ C^T
 * K C \f$ and \f$ C^T F \f$ on their unconstrained rows and columns, so dropping their constrained
 * rows and columns at insertion assembles the reduced system exactly.
 */
class ConstraintEliminationDofMap : public ReducedDofMap
{
public:
  ConstraintEliminationDofMap(MeshBase & mesh, System & system, const DofMap & dof_map);
  virtual ~ConstraintEliminationDofMap();

  /**
   * Unconstrained degrees of freedom are retained
   */
  virtual bool retain_field_dof(const Elem & elem,
                                unsigned int node_num,
                                unsigned int var_num,
                                dof_id_type field_dof) const override;

  /**
   * Unconstrained degrees of freedom are retained
   */
  virtual bool retain_scalar_dof(dof_id_type scalar_dof) const override;

  /**
   * Prolongs \p reduced, a vector on this space, into the full-space vector \p full, which must be
   * \p PARALLEL or \p GHOSTED. The ghost entries of a \p GHOSTED vector must include the
   * \p DofMap send list, and they are updated. With \p homogeneous false this forms
   * \f$ x = C y + h \f$; with \p homogeneous true, as for a step, a search direction, or an
   * eigenvector, \f$ x = C y \f$
   */
  void prolong(const NumericVector<Number> & reduced,
               NumericVector<Number> & full,
               bool homogeneous) const;

private:
  bool is_unconstrained(dof_id_type dof) const;
};

}

#endif // LIBMESH_CONSTRAINT_ELIMINATION_DOF_MAP_H
