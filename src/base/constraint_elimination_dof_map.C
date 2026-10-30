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

#include "libmesh/constraint_elimination_dof_map.h"

#include "libmesh/dof_map.h"

namespace libMesh
{
ConstraintEliminationDofMap::ConstraintEliminationDofMap(MeshBase & mesh,
                                                         System & system,
                                                         const DofMap & dof_map)
  : ReducedDofMap(mesh, system, dof_map, "constraint_eliminated")
{
  libmesh_experimental();
}

ConstraintEliminationDofMap::~ConstraintEliminationDofMap() = default;

bool ConstraintEliminationDofMap::retain_field_dof(const Elem &,
                                                   unsigned int,
                                                   unsigned int,
                                                   const dof_id_type field_dof) const
{
  return this->is_unconstrained(field_dof);
}

bool ConstraintEliminationDofMap::retain_scalar_dof(const dof_id_type scalar_dof) const
{
  return this->is_unconstrained(scalar_dof);
}

bool ConstraintEliminationDofMap::is_unconstrained(const dof_id_type dof) const
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  return !_dof_map.is_constrained_dof(dof);
#else
  libmesh_ignore(dof);
  return true;
#endif
}
}
