// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_dm_wrapper.h"

#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

namespace libMesh
{

void PetscDMWrapper::set_point_range_in_section( const System & system, PetscSection & section)
{
  const MeshBase & mesh = system.get_mesh();

  unsigned int pstart = 2^30;
  unsigned int pend = 0;

  // Find minimum and maximum (inclusive) global id numbers on the current processor
  // to build PetscSection. Currently, we're using the global node number as the point index.
  // TODO: This is currently restricted to nodal dofs only!
  //       Need to generalize to the cases where there's dofs at edges, faces, and interiors
  //       as well as vertices. When we generalize, we must guarantee that each id() is unique!
  //       We'll then add each edge, face, and/or interior as a new "point" in the PetscSection.
  if (mesh.n_active_local_elem() > 0)
    {
      for (MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
           el != mesh.active_local_elements_end(); el++)
        {
          const Elem * elem = *el;

          if (mesh.query_elem_ptr(elem->id()))
            {
              for (unsigned int n = 0; n < elem->n_nodes(); n++)
                {
                  // get the global id number of local node n
                  dof_id_type node = elem->node(n);

                  if (node < pstart)
                    pstart = node;

                  if (node > pend)
                    pend = node;
                }
            }
        }

      // PetscSectionSetChart is expecting [pstart,pend), so pad pend by 1.
      pend +=1;
    }

  // If we're on a processor who coarsened the mesh to have no local elements,
  // we should make an empty PetscSection. An empty PetscSection is specified
  // by passing [0,0) to the PetscSectionSetChart call.
  else
    {
      pstart = 0;
      pend = 0;
    }

  PetscErrorCode ierr = PetscSectionSetChart(section, pstart, pend);
  CHKERRABORT(system.comm().get(),ierr);
}

} // end namespace libMesh

#endif // LIBMESH_HAVE_PETSC
