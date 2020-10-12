// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/ignore_warnings.h"
#include <petscsf.h>
#include "libmesh/restore_warnings.h"

#include "libmesh/petsc_dm_wrapper.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"

namespace libMesh
{

PetscDMWrapper::~PetscDMWrapper()
{
  this->clear();
}

void PetscDMWrapper::clear()
{
  // This will also Destroy the attached PetscSection and PetscSF as well
  // Destroy doesn't free the memory, but just resets points internally
  // in the struct, so we'd still need to wipe out the memory on our side
  for( auto dm_it = _dms.begin(); dm_it < _dms.end(); ++dm_it )
    DMDestroy( dm_it->get() );

  _dms.clear();
  _sections.clear();
  _star_forests.clear();
}

void PetscDMWrapper::init_and_attach_petscdm(System & system, SNES & snes)
{
  START_LOG ("init_and_attach_petscdm()", "PetscDMWrapper");

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
#if PETSC_VERSION_LESS_THAN(3,10,0)
      ierr = DMSetDefaultSection(dm, section);
#else
      ierr = DMSetSection(dm, section);
#endif
      CHKERRABORT(system.comm().get(),ierr);

      // We only need to build the star forest if we're in a parallel environment
      if (system.n_processors() > 1)
        {
          // Build the PetscSF and attach it to the DM
          this->build_sf(system, star_forest);
#if PETSC_VERSION_LESS_THAN(3,12,0)
          ierr = DMSetDefaultSF(dm, star_forest);
#else
          ierr = DMSetSectionSF(dm, star_forest);
#endif
          CHKERRABORT(system.comm().get(),ierr);
        }
    }

  // We need to set only the finest level DM in the SNES
  DM & dm = this->get_dm(0);
  ierr = SNESSetDM(snes, dm);
  CHKERRABORT(system.comm().get(),ierr);

  STOP_LOG ("init_and_attach_petscdm()", "PetscDMWrapper");
}

void PetscDMWrapper::build_section( const System & system, PetscSection & section )
{
  START_LOG ("build_section()", "PetscDMWrapper");

  PetscErrorCode ierr;
  ierr = PetscSectionCreate(system.comm().get(),&section);
  CHKERRABORT(system.comm().get(),ierr);

  // Tell the PetscSection about all of our System variables
  ierr = PetscSectionSetNumFields(section,system.n_vars());
  CHKERRABORT(system.comm().get(),ierr);

  // Set the actual names of all the field variables
  for( unsigned int v = 0; v < system.n_vars(); v++ )
    {
      ierr = PetscSectionSetFieldName( section, v, system.variable_name(v).c_str() );
      CHKERRABORT(system.comm().get(),ierr);
    }

  // For building the section, we need to create local-to-global map
  // of local "point" ids to the libMesh global id of that point.
  // A "point" in PETSc nomenclature is a geometric object that can have
  // dofs associated with it, e.g. Node, Edge, Face, Elem.
  // The numbering PETSc expects is continuous for the local numbering.
  // Since we're only using this interface for solvers, then we can just
  // assign whatever local id to any of the global ids. But it is local
  // so we don't need to worry about processor numbering for the local
  // point ids.
  std::unordered_map<dof_id_type,dof_id_type> node_map;
  std::unordered_map<dof_id_type,dof_id_type> elem_map;
  std::map<dof_id_type,unsigned int> scalar_map;

  // First we tell the PetscSection about all of our points that have
  // dofs associated with them.
  this->set_point_range_in_section(system, section, node_map, elem_map, scalar_map);

  // Now we can build up the dofs per "point" in the PetscSection
  this->add_dofs_to_section(system, section, node_map, elem_map, scalar_map);

  // Final setup of PetscSection
  // Until Matt Knepley finishes implementing the commented out function
  // below, the PetscSection will be assuming node-major ordering
  // so let's throw and error if the user tries to use this without
  // node-major order
  if (!libMesh::on_command_line("--node-major-dofs"))
    libmesh_error_msg("ERROR: Must use --node-major-dofs with PetscSection!");

  //else if (!system.identify_variable_groups())
  //  ierr = PetscSectionSetUseFieldOffsets(section,PETSC_TRUE);CHKERRABORT(system.comm().get(),ierr);
  //else
  //  {
  //    std::string msg = "ERROR: Only node-major or var-major ordering supported for PetscSection!\n";
  //    msg += "       var-group-major ordering not supported!\n";
  //    msg += "       Must use --node-major-dofs or set System::identify_variable_groups() = false!\n";
  //    libmesh_error_msg(msg);
  //  }

  ierr = PetscSectionSetUp(section);CHKERRABORT(system.comm().get(),ierr);

  // Sanity checking at least that local_n_dofs match
  libmesh_assert_equal_to(system.n_local_dofs(),this->check_section_n_dofs(system, section));

  STOP_LOG ("build_section()", "PetscDMWrapper");
}

void PetscDMWrapper::build_sf( const System & system, PetscSF & star_forest )
{
  START_LOG ("build_sf()", "PetscDMWrapper");

  const DofMap & dof_map = system.get_dof_map();

  const std::vector<dof_id_type> & send_list = dof_map.get_send_list();

  // Number of ghost dofs that send information to this processor
  const PetscInt n_leaves = cast_int<PetscInt>(send_list.size());

  // Number of local dofs, including ghosts dofs
  const PetscInt n_roots = dof_map.n_local_dofs() + n_leaves;

  // This is the vector of dof indices coming from other processors
  // We need to give this to the PetscSF
  // We'll be extra paranoid about this ugly double cast
  static_assert(sizeof(PetscInt) == sizeof(dof_id_type),"PetscInt is not a dof_id_type!");
  PetscInt * local_dofs = reinterpret_cast<PetscInt *>(const_cast<dof_id_type *>(send_list.data()));

  // This is the vector of PetscSFNode's for the local_dofs.
  // For each entry in local_dof, we have to supply the rank from which
  // that dof stems and its local index on that rank.
  // PETSc documentation here:
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PetscSF/PetscSFNode.html
  std::vector<PetscSFNode> sf_nodes(send_list.size());

  for( unsigned int i = 0; i < send_list.size(); i++ )
    {
      dof_id_type incoming_dof = send_list[i];

      const processor_id_type rank = dof_map.dof_owner(incoming_dof);

      // Dofs are sorted and continuous on the processor so local index
      // is counted up from the first dof on the processor.
      PetscInt index = incoming_dof - dof_map.first_dof(rank);

      sf_nodes[i].rank  = rank; /* Rank of owner */
      sf_nodes[i].index = index;/* Index of dof on rank */
    }

  PetscSFNode * remote_dofs = sf_nodes.data();

  PetscErrorCode ierr;
  ierr = PetscSFCreate(system.comm().get(), &star_forest);CHKERRABORT(system.comm().get(),ierr);

  // TODO: We should create pointers to arrays so we don't have to copy
  //       and then can use PETSC_OWN_POINTER where PETSc will take ownership
  //       and delete the memory for us. But then we'd have to use PetscMalloc.
  ierr = PetscSFSetGraph(star_forest,
                         n_roots,
                         n_leaves,
                         local_dofs,
                         PETSC_COPY_VALUES,
                         remote_dofs,
                         PETSC_COPY_VALUES);
  CHKERRABORT(system.comm().get(),ierr);

  STOP_LOG ("build_sf()", "PetscDMWrapper");
}

void PetscDMWrapper::set_point_range_in_section (const System & system,
                                                 PetscSection & section,
                                                 std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                                 std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                                                 std::map<dof_id_type,unsigned int> & scalar_map)
{
  // We're expecting this to be empty coming in
  libmesh_assert(node_map.empty());

  // We need to count up the number of active "points" on this processor.
  // Nominally, a "point" in PETSc parlance is a geometric object that can
  // hold DoFs, i.e node, edge, face, elem. Since we handle the mesh and are only
  // interested in solvers, then the only thing PETSc needs is a unique *local* number
  // for each "point" that has active DoFs; note however this local numbering
  // we construct must be continuous.
  //
  // In libMesh, for most finite elements, we just associate those DoFs with the
  // geometric nodes. So can we loop over the nodes on this processor and check
  // if any of the fields are have active DoFs on that node.
  // If so, then we tell PETSc about that "point". At this stage, we just need
  // to count up how many active "points" we have and cache the local number to global id
  // mapping.


  // These will be our local counters. pstart should always be zero.
  // pend will track our local "point" count.
  // If we're on a processor who coarsened the mesh to have no local elements,
  // we should make an empty PetscSection. An empty PetscSection is specified
  // by passing [0,0) to the PetscSectionSetChart call at the end. So, if we
  // have nothing on this processor, these are the correct values to pass to
  // PETSc.
  dof_id_type pstart = 0;
  dof_id_type pend = 0;

  const MeshBase & mesh = system.get_mesh();

  const DofMap & dof_map = system.get_dof_map();

  // If we don't have any local dofs, then there's nothing to tell to the PetscSection
  if (dof_map.n_local_dofs() > 0)
    {
      // Conservative estimate of space needed so we don't thrash
      node_map.reserve(mesh.n_local_nodes());
      elem_map.reserve(mesh.n_active_local_elem());

      // We loop over active elements and then cache the global/local node mapping to make sure
      // we only count active nodes. For example, if we're calling this function and we're
      // not the finest level in the Mesh tree, we don't want to include nodes of child elements
      // that aren't active on this level.
      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
            {
              // get the global id number of local node n
              const Node & node = elem->node_ref(n);

              // Only register nodes with the PetscSection if they have dofs that belong to
              // this processor. Even though we're active local elements, the dofs associated
              // with the node may belong to a different processor. The processor who owns
              // those dofs will register that node with the PetscSection on that processor.
              std::vector<dof_id_type> node_dof_indices;
              dof_map.dof_indices( &node, node_dof_indices );
              if( !node_dof_indices.empty() && dof_map.local_index(node_dof_indices[0]) )
                {
#ifndef NDEBUG
                  // We're assuming that if the first dof for this node belongs to this processor,
                  // then all of them do.
                  for( auto dof : node_dof_indices )
                    libmesh_assert(dof_map.local_index(dof));
#endif
                  // Cache the global/local mapping if we haven't already
                  // Then increment our local count
                  dof_id_type node_id = node.id();
                  if( node_map.count(node_id) == 0 )
                    {
                      node_map.insert(std::make_pair(node_id,pend));
                      pend++;
                    }
                }
            }

          // Some finite elements, e.g. Hierarchic, associate element interior DoFs with the element
          // rather than the node (since we ought to be able to use Hierachic elements on a QUAD4,
          // which has no interior node). Thus, we also need to check element interiors for DoFs
          // as well and, if the finite element has them, we also need to count the Elem in our
          // "point" accounting.
          if( elem->n_dofs(system.number()) > 0 )
            {
              dof_id_type elem_id = elem->id();
              elem_map.insert(std::make_pair(elem_id,pend));
              pend++;
            }
        }

      // SCALAR dofs live on the "last" processor, so only work there if there are any
      if( dof_map.n_SCALAR_dofs() > 0 && (system.processor_id() == (system.n_processors()-1)) )
        {
          // Loop through all the variables and cache the scalar ones. We cache the
          // SCALAR variable index along with the local point to make it easier when
          // we have to register dofs with the PetscSection
          for( unsigned int v = 0; v < system.n_vars(); v++ )
            {
              if( system.variable(v).type().family == SCALAR )
                {
                  scalar_map.insert(std::make_pair(pend,v));
                  pend++;
                }
            }
        }

    }

  PetscErrorCode ierr = PetscSectionSetChart(section, pstart, pend);
  CHKERRABORT(system.comm().get(),ierr);
}

void PetscDMWrapper::add_dofs_to_section (const System & system,
                                          PetscSection & section,
                                          const std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                          const std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                                          const std::map<dof_id_type,unsigned int> & scalar_map)
{
  const MeshBase & mesh = system.get_mesh();

  PetscErrorCode ierr;

  // Now we go through and add dof information for each "point".
  //
  // In libMesh, for most finite elements, we just associate those DoFs with the
  // geometric nodes. So can we loop over the nodes we cached in the node_map and
  // the DoFs for each field for that node. We need to give PETSc the local id
  // we built up in the node map.
  for (const auto & nmap : node_map )
    {
      const dof_id_type global_node_id = nmap.first;
      const dof_id_type local_node_id = nmap.second;

      libmesh_assert( mesh.query_node_ptr(global_node_id) );

      const Node & node = mesh.node_ref(global_node_id);

      this->add_dofs_helper(system,node,local_node_id,section);
    }

  // Some finite element types associate dofs with the element. So now we go through
  // any of those with the Elem as the point we add to the PetscSection with accompanying
  // dofs
  for (const auto & emap : elem_map )
    {
      const dof_id_type global_elem_id = emap.first;
      const dof_id_type local_elem_id = emap.second;

      libmesh_assert( mesh.query_elem_ptr(global_elem_id) );

      const Elem & elem = mesh.elem_ref(global_elem_id);

      this->add_dofs_helper(system,elem,local_elem_id,section);
    }

  // Now add any SCALAR dofs to the PetscSection
  // SCALAR dofs live on the "last" processor, so only work there if there are any
  if (system.processor_id() == (system.n_processors()-1))
    {
      for (const auto & smap : scalar_map )
        {
          const dof_id_type local_id = smap.first;
          const unsigned int scalar_var = smap.second;

          // The number of SCALAR dofs comes from the variable order
          const int n_dofs = system.variable(scalar_var).type().order.get_order();

          ierr = PetscSectionSetFieldDof( section, local_id, scalar_var, n_dofs );
          CHKERRABORT(system.comm().get(),ierr);

          // In the SCALAR case, there are no other variables associate with the "point"
          // the total number of dofs on the point is the same as that for the field
          ierr = PetscSectionSetDof( section, local_id, n_dofs );
          CHKERRABORT(system.comm().get(),ierr);
        }
    }

}

void PetscDMWrapper::add_dofs_helper (const System & system,
                                      const DofObject & dof_object,
                                      dof_id_type local_id,
                                      PetscSection & section)
{
  unsigned int total_n_dofs_at_dofobject = 0;

  // We are assuming variables are also numbered 0 to n_vars()-1
  for( unsigned int v = 0; v < system.n_vars(); v++ )
    {
      unsigned int n_dofs_at_dofobject = dof_object.n_dofs(system.number(), v);

      if( n_dofs_at_dofobject > 0 )
        {
          PetscErrorCode ierr = PetscSectionSetFieldDof( section,
                                                         local_id,
                                                         v,
                                                         n_dofs_at_dofobject );

          CHKERRABORT(system.comm().get(),ierr);

          total_n_dofs_at_dofobject += n_dofs_at_dofobject;
        }
    }

  libmesh_assert_equal_to(total_n_dofs_at_dofobject, dof_object.n_dofs(system.number()));

  PetscErrorCode ierr =
    PetscSectionSetDof( section, local_id, total_n_dofs_at_dofobject );
  CHKERRABORT(system.comm().get(),ierr);
}


dof_id_type PetscDMWrapper::check_section_n_dofs( const System & system, PetscSection & section )
{
  PetscInt n_local_dofs = 0;

  // Grap the starting and ending points from the section
  PetscInt pstart, pend;
  PetscErrorCode ierr = PetscSectionGetChart(section, &pstart, &pend);
  CHKERRABORT(system.comm().get(),ierr);

  // Count up the n_dofs for each point from the section
  for( PetscInt p = pstart; p < pend; p++ )
    {
      PetscInt n_dofs;
      ierr = PetscSectionGetDof(section,p,&n_dofs);CHKERRABORT(system.comm().get(),ierr);
      n_local_dofs += n_dofs;
    }

  static_assert(sizeof(PetscInt) == sizeof(dof_id_type),"PetscInt is not a dof_id_type!");
  return n_local_dofs;
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
