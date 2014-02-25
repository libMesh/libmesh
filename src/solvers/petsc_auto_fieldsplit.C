// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#if !PETSC_VERSION_LESS_THAN(3,2,0)

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/system.h"

EXTERN_C_FOR_PETSC_BEGIN
#  include <petscksp.h>
EXTERN_C_FOR_PETSC_END

// C++ includes

namespace libMesh
{

void petsc_auto_fieldsplit (PC my_pc,
                            const System &sys)
{
  if (libMesh::on_command_line("--solver_variable_names"))
    {
      for (unsigned int v = 0; v != sys.n_vars(); ++v)
        {
          const std::string& var_name = sys.variable_name(v);
          std::vector<dof_id_type> var_idx;
          sys.get_dof_map().local_variable_indices
            (var_idx, sys.get_mesh(), v);

          IS is;

          PetscInt *idx = PETSC_NULL;
          if (!var_idx.empty())
            idx = reinterpret_cast<PetscInt*>(&var_idx[0]);

          int ierr = ISCreateLibMesh(sys.comm().get(), var_idx.size(),
                                     idx, PETSC_COPY_VALUES, &is);
          CHKERRABORT(sys.comm().get(), ierr);

          ierr = PCFieldSplitSetIS(my_pc, var_name.c_str(), is);
          CHKERRABORT(sys.comm().get(), ierr);
        }
    }
}

} // namespace libMesh


#else  // #PETSC_VERSION < 3.2.0
void assign_solver_fieldsplit_names (PC my_pc,
                                     const System &sys)
{
  if (libMesh::on_command_line("--solver_variable_names"))
    {
      libmesh_do_once(
        libMesh::out <<
          "WARNING: libMesh does not support setting field splits" <<
          std::endl << "with PETSc " <<
          LIBMESH_DETECTED_PETSC_VERSION_MAJOR << '.'
          LIBMESH_DETECTED_PETSC_VERSION_MINOR << '.'
          LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR << std::endl;);
    }
}
#endif // #PETSC_VERSION > 3.2.0
#endif // #ifdef LIBMESH_HAVE_PETSC
