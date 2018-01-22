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

// PETSc includes
#include <petscsf.h>

#include "libmesh/petsc_dm_wrapper.h"

#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"

namespace libMesh
{

void PetscDMWrapper::init_and_attach_petscdm(const System & system, SNES & snes)
{
  START_LOG ("init_and_attach_petscdm", "PetscDMWrapper");

  PetscErrorCode ierr;

  // Eventually, we'll traverse the mesh hierarchy and cache
  // for each grid level, but for now, we're just using the
  // finest level
  unsigned int n_levels = 1;
  this->init_dm_data(n_levels);

  for(unsigned int level = 0; level < n_levels; level++)
    {
      DM & dm = this->get_dm(level);
      PetscSection & section = this->get_section(level);
      PetscSF & star_forest = this->get_star_forest(level);

      ierr = DMShellCreate(system.comm().get(), &dm);
      CHKERRABORT(system.comm().get(),ierr);

      // Build the PetscSection and attach it to the DM
      this->build_section(system, section);
      ierr = DMSetDefaultSection(dm, section);
      CHKERRABORT(system.comm().get(),ierr);

      // We only need to build the star forest if we're in a parallel environment
      if (system.n_processors() > 1)
        {
          // Build the PetscSF and attach it to the DM
          this->build_sf(system, star_forest);
          ierr = DMSetDefaultSF(dm, star_forest);
          CHKERRABORT(system.comm().get(),ierr);
        }
    }

  // We need to set only the finest level DM in the SNES
  DM & dm = this->get_dm(0);
  ierr = SNESSetDM(snes, dm);
  CHKERRABORT(system.comm().get(),ierr);

  STOP_LOG ("init_and_attach_petscdm", "PetscDMWrapper");
}

void PetscDMWrapper::build_section( const System & system, PetscSection & section )
{
  START_LOG ("build_section", "PetscDMWrapper");

  PetscErrorCode ierr;
  ierr = PetscSectionCreate(system.comm().get(),&section);
  CHKERRABORT(system.comm().get(),ierr);

  ierr = PetscSectionSetNumFields(section,system.n_vars());
  CHKERRABORT(system.comm().get(),ierr);

  // First, set the actual names of all the fields variables we are interested in
  for( unsigned int v = 0; v < system.n_vars(); v++ )
    {
      ierr = PetscSectionSetFieldName( section, v, system.variable_name(v).c_str() );
      CHKERRABORT(system.comm().get(),ierr);
    }

  // Set "points" count into the section. A "point" in PETSc nomenclature
  // is a geometric object that can have dofs associated with it, e.g.
  // Node, Edge, Face, Elem. First we tell the PetscSection about all of our
  // points.
  this->set_point_range_in_section(system, section);

  // Now build up the dofs per "point" in the PetscSection
  this->add_dofs_to_section(system, section);

  // Final setup of PetscSection
  ierr = PetscSectionSetUp(section);CHKERRABORT(system.comm().get(),ierr);

  // Sanity checking at least that total n_dofs match
#ifndef NDEBUG
  this->check_section_n_dofs_match(system, section);
#endif

  STOP_LOG ("build_section", "PetscDMWrapper");
}

void PetscDMWrapper::build_sf( const System & system, PetscSF & star_forest )
{
  START_LOG ("build_sf", "PetscDMWrapper");

  const DofMap & dof_map = system.get_dof_map();

  const std::vector<dof_id_type> & send_list = dof_map.get_send_list();

  // Number of ghost dofs that send information to this processor
  PetscInt n_leaves = send_list.size();

  // Number of local dofs, including ghosts dofs
  PetscInt n_roots = dof_map.n_local_dofs();
  n_roots += send_list.size();

  // This is the vector of dof indices coming from other processors
  // TODO: We do a stupid copy here since we can't convert an
  //       unsigned int* to PetscInt*
  std::vector<PetscInt> send_list_copy(send_list.size());
  for( unsigned int i = 0; i < send_list.size(); i++ )
    send_list_copy[i] = send_list[i];

  PetscInt * local_dofs = send_list_copy.data();

  // This is the vector of PetscSFNode's for the local_dofs.
  // For each entry in local_dof, we have to supply the rank from which
  // that dof stems and its local index on that rank.
  // PETSc documentation here:
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PetscSF/PetscSFNode.html
  std::vector<PetscSFNode> sf_nodes(send_list.size());

  for( unsigned int i = 0; i < send_list.size(); i++ )
    {
      unsigned int incoming_dof = send_list[i];

      PetscInt rank = this->find_dof_rank( incoming_dof, dof_map );

      /* Dofs are sorted and continuous on the processor so local index
         is counted up from the first dof on the processor. */
      PetscInt index = incoming_dof - dof_map.first_dof(rank);

      sf_nodes[i].rank  = rank; /* Rank of owner */
      sf_nodes[i].index = index;/* Index of dof on rank */
    }

  PetscSFNode * remote_dofs = sf_nodes.data();

  PetscErrorCode ierr;
  ierr = PetscSFCreate(system.comm().get(), &star_forest);CHKERRABORT(system.comm().get(),ierr);

  // TODO: We should create pointers to arrays so we don't have to copy
  //       and case use PETSC_OWN_POINTER where PETSc will take ownership
  //       and delete the memory for us. But then we have to use PetscMalloc.
  ierr = PetscSFSetGraph(star_forest,
                         n_roots,
                         n_leaves,
                         local_dofs,
                         PETSC_COPY_VALUES,
                         remote_dofs,
                         PETSC_COPY_VALUES);
  CHKERRABORT(system.comm().get(),ierr);

  STOP_LOG ("build_sf", "PetscDMWrapper");
}

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

void PetscDMWrapper::add_dofs_to_section( const System & system, PetscSection & section )
{
  const MeshBase & mesh = system.get_mesh();
  const DofMap & dof_map = system.get_dof_map();

  PetscErrorCode ierr;

  // Now we go through and add dof information for each point object.
  // We do this by looping through all the elements on this processor
  // and add associated dofs with that element
  // TODO: Currently, we're only adding the dofs at nodes! We need to generalize
  //       beyond nodal FEM!
  for (MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
       el != mesh.active_local_elements_end(); el++)
    {
      const Elem * elem = *el;

      if (mesh.query_elem_ptr(elem->id()))
        {
          // We need to keep a count of total dofs at each point
          // for the PetscSectionSetDof call at the end.
          std::map<unsigned int, unsigned int> global_dofs;

          // Now set the dof for each field
          for( unsigned int v = 0; v < system.n_vars(); v++ )
            {
              std::vector<dof_id_type> dof_indices;
              dof_map.dof_indices( elem, dof_indices, v );

              // Right now, we're assuming at most one dof per node
              // TODO:  Need to remove this assumption of at most one dof per node
              for( unsigned int n = 0; n < dof_indices.size(); n++ )
                {
                  dof_id_type index = dof_indices[n];

                  // Check if this dof index is on this processor. If so, we count it.
                  if( index >= dof_map.first_dof() && index < dof_map.end_dof() )
                    {
                      // TODO: Need to remove this assumption of at most one dof per node
                      PetscInt n_dofs_at_node = 1;

                      dof_id_type global_node = elem->node(n);

                      ierr = PetscSectionSetFieldDof( section, global_node, v, n_dofs_at_node );
                      CHKERRABORT(system.comm().get(),ierr);

                      global_dofs[global_node] += 1;
                    }
                }
            }

          // [PB]: This is redundant, but PETSc needed it at the time of writing. Perhaps
          //       it will be fixed upstream at some point.
          for( std::map<unsigned int, unsigned int>::const_iterator it = global_dofs.begin();
               it != global_dofs.end(); ++it )
            {
              unsigned int global_node = it->first;
              unsigned int total_n_dofs_at_node = it->second;

              ierr = PetscSectionSetDof( section, global_node, total_n_dofs_at_node );
              CHKERRABORT(system.comm().get(),ierr);
            }
        }
    }
}

void PetscDMWrapper::check_section_n_dofs_match( const System & system, PetscSection & section )
{
  PetscInt total_n_dofs = 0;

  // Grap the starting and ending points from the section
  PetscInt pstart, pend;
  PetscErrorCode ierr = PetscSectionGetChart(section, &pstart, &pend);
  CHKERRABORT(system.comm().get(),ierr);

  // Count up the n_dofs for each point from the section
  for( PetscInt p = pstart; p < pend+1; p++ )
    {
      PetscInt n_dofs;
      ierr = PetscSectionGetDof(section,p,&n_dofs);CHKERRABORT(system.comm().get(),ierr);
      total_n_dofs += n_dofs;
    }

  // That should match the n_local_dofs for our system
  libmesh_assert_equal_to(total_n_dofs,(PetscInt)system.n_local_dofs());
}

PetscInt PetscDMWrapper::find_dof_rank( unsigned int dof, const DofMap& dof_map ) const
{
  libmesh_assert_greater_equal( dof, dof_map.first_dof( 0 ) );
  libmesh_assert_less( dof, dof_map.end_dof(dof_map.comm().size()-1) );

  // dofs are in order on each processor starting from processor 0, so we can
  // use a binary search
  unsigned int max_rank = dof_map.comm().size();
  unsigned int current_rank = max_rank/2;
  unsigned int min_rank = 0;

  do
    {
      if( dof >= dof_map.first_dof(current_rank) )
        {
          min_rank = current_rank;
          current_rank = (max_rank + min_rank)/2;
        }
      else
        {
          max_rank = current_rank;
          current_rank = (max_rank + min_rank)/2;
        }
    }
  while( max_rank - min_rank > 1 );

  libmesh_assert_less( dof, dof_map.end_dof(current_rank) );
  libmesh_assert_greater_equal( dof, dof_map.first_dof(current_rank) );

  return current_rank;
}

void PetscDMWrapper::init_dm_data(unsigned int n_levels)
{
  _dms.resize(n_levels);
  _sections.resize(n_levels);
  _star_forests.resize(n_levels);

  for( unsigned int i = 0; i < n_levels; i++ )
    {
      _dms[i] = libmesh_make_unique<DM>();
      _sections[i] = libmesh_make_unique<PetscSection>();
      _star_forests[i] = libmesh_make_unique<PetscSF>();
    }
}

} // end namespace libMesh

#endif // LIBMESH_HAVE_PETSC
