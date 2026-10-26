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
#include "libmesh/numeric_vector.h"
#include "libmesh/int_range.h"

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

void ConstraintEliminationDofMap::prolong(const NumericVector<Number> & reduced,
                                          NumericVector<Number> & full,
                                          const bool homogeneous) const
{
  // Each rank writes the retained entries it owns. A serial vector holds every entry on every rank,
  // so the entries owned by other ranks would keep stale values
  libmesh_assert_not_equal_to(full.type(), SERIAL);
  libmesh_assert_equal_to(reduced.local_size(), this->n_local_dofs());

  std::vector<numeric_index_type> reduced_indices(this->n_local_dofs());
  for (const auto i : index_range(reduced_indices))
    reduced_indices[i] = this->first_dof() + i;
  std::vector<Number> values;
  reduced.get(reduced_indices, values);

  // The constrained entries of full are overwritten by the constraints, so their prior values do
  // not matter
  full.insert(values, this->local_retained_dofs());
  full.close();
  _dof_map.enforce_constraints_exactly(_system, &full, homogeneous);
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
