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
#include <vector>

namespace libMesh
{
class Variable;
class Elem;
class Node;

class DofMapBase
{
public:
  DofMapBase() = default;

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
  virtual dof_id_type first_dof() const = 0;

  /**
   * \returns The first dof index that is after all indices local to
   * processor \p proc.
   *
   * Analogous to the end() member function of STL containers.
   */
  virtual dof_id_type end_dof() const = 0;

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
  virtual dof_id_type n_dofs() const = 0;

  /**
   * \returns The number of degrees of freedom on this processor.
   */
  virtual dof_id_type n_local_dofs() const = 0;
};

}
#endif // LIBMESH_DOF_MAP_BASE_H
