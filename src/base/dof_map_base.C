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

#include "libmesh/dof_map_base.h"
#include "libmesh/parallel_implementation.h"

namespace libMesh
{
DofMapBase::DofMapBase(const Parallel::Communicator & comm)
  : ParallelObject(comm),
    _n_dfs(0)
#ifdef LIBMESH_ENABLE_AMR
    ,
    _n_old_dfs(0)
#endif
{
}

std::size_t DofMapBase::compute_dof_info(const dof_id_type n_local_dofs)
{
  // Get DOF counts on all processors
  const auto n_proc = this->n_processors();

  std::vector<dof_id_type> dofs_on_proc(n_proc, 0);
  this->comm().allgather(n_local_dofs, dofs_on_proc);

#ifdef LIBMESH_ENABLE_AMR
  // Resize and fill the _first_df and _end_df arrays
  _first_old_df = _first_df;
  _end_old_df = _end_df;
#endif

  _first_df.resize(n_proc);
  _end_df.resize(n_proc);

  // Get DOF offsets
  _first_df[0] = 0;
  for (processor_id_type i = 1; i < n_proc; ++i)
    _first_df[i] = _end_df[i - 1] = _first_df[i - 1] + dofs_on_proc[i - 1];
  _end_df[n_proc - 1] = _first_df[n_proc - 1] + dofs_on_proc[n_proc - 1];

// Set the total number of degrees of freedom
#ifdef LIBMESH_ENABLE_AMR
  _n_old_dfs = _n_dfs;
#endif
  _n_dfs = _end_df[n_proc - 1];

  // Return total number of DOFs across all procs. We compute and
  // return this as a std::size_t so that we can detect situations in
  // which the total number of DOFs across all procs would exceed the
  // capability of the underlying NumericVector representation to
  // index into it correctly (std::size_t is the largest unsigned
  // type, so no NumericVector representation can exceed it).
  return std::accumulate(dofs_on_proc.begin(), dofs_on_proc.end(), static_cast<std::size_t>(0));
}

void DofMapBase::clear()
{
  _first_df.clear();
  _end_df.clear();
  _n_dfs = 0;
}
}
