// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/petsc_auto_fieldsplit.h"

#ifdef LIBMESH_HAVE_PETSC

#include <petscksp.h>

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/system.h"

namespace {
using namespace libMesh;

void indices_to_fieldsplit (const Parallel::Communicator & comm,
                            const std::vector<dof_id_type> & indices,
                            PC my_pc,
                            const std::string & field_name)
{
  const PetscInt * idx = PETSC_NULL;
  if (!indices.empty())
    idx = reinterpret_cast<const PetscInt *>(indices.data());

  IS is;
  int ierr = ISCreateGeneral(comm.get(), cast_int<PetscInt>(indices.size()),
                             idx, PETSC_COPY_VALUES, &is);
  CHKERRABORT(comm.get(), ierr);

  ierr = PCFieldSplitSetIS(my_pc, field_name.c_str(), is);
  CHKERRABORT(comm.get(), ierr);
}

}

namespace libMesh
{

void petsc_auto_fieldsplit (PC my_pc,
                            const System & sys)
{
  std::string sys_prefix = "--solver_group_";

  if (libMesh::on_command_line("--solver-system-names"))
    {
      sys_prefix = sys_prefix + sys.name() + "_";
    }

  std::map<std::string, std::vector<dof_id_type>> group_indices;

  if (libMesh::on_command_line("--solver-variable-names"))
    {
      for (auto v : make_range(sys.n_vars()))
        {
          const std::string & var_name = sys.variable_name(v);

          std::vector<dof_id_type> var_idx;
          sys.get_dof_map().local_variable_indices
            (var_idx, sys.get_mesh(), v);

          std::string group_command = sys_prefix + var_name;

          const std::string empty_string;

          std::string group_name = libMesh::command_line_value
            (group_command, empty_string);

          if (group_name != empty_string)
            {
              std::vector<dof_id_type> & indices =
                group_indices[group_name];
              const bool prior_indices = !indices.empty();
              indices.insert(indices.end(), var_idx.begin(),
                             var_idx.end());
              if (prior_indices)
                std::sort(indices.begin(), indices.end());
            }
          else
            {
              indices_to_fieldsplit (sys.comm(), var_idx, my_pc, var_name);
            }
        }
    }

  for (const auto & pr : group_indices)
    indices_to_fieldsplit(sys.comm(), pr.second, my_pc, pr.first);
}

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
