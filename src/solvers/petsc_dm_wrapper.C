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

#include <petscsf.h>

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

  // First we tell the PetscSection about all of our points that have
  // dofs associated with them.
  this->set_point_range_in_section(system, section, node_map, elem_map);

  // Now we can build up the dofs per "point" in the PetscSection
  this->add_dofs_to_section(system, section, node_map, elem_map);

  // Final setup of PetscSection
  ierr = PetscSectionSetUp(section);CHKERRABORT(system.comm().get(),ierr);

  STOP_LOG ("build_section()", "PetscDMWrapper");
}

void PetscDMWrapper::set_point_range_in_section (const System & system,
                                                 PetscSection & section,
                                                 std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                                 std::unordered_map<dof_id_type,dof_id_type> & elem_map)
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

    }

  PetscErrorCode ierr = PetscSectionSetChart(section, pstart, pend);
  CHKERRABORT(system.comm().get(),ierr);
}

void PetscDMWrapper::add_dofs_to_section (const System & system,
                                          PetscSection & section,
                                          const std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                          const std::unordered_map<dof_id_type,dof_id_type> & elem_map)
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

      unsigned int total_n_dofs_at_node = 0;

      // We are assuming variables are also numbered 0 to n_vars()-1
      for( unsigned int v = 0; v < system.n_vars(); v++ )
        {
          unsigned int n_dofs_at_node = node.n_dofs(system.number(), v);

          if( n_dofs_at_node > 0 )
            {
              ierr = PetscSectionSetFieldDof( section, local_node_id, v, n_dofs_at_node );
              CHKERRABORT(system.comm().get(),ierr);

              total_n_dofs_at_node += n_dofs_at_node;
            }
        }

      libmesh_assert_equal_to(total_n_dofs_at_node, node.n_dofs(system.number()));

      ierr = PetscSectionSetDof( section, local_node_id, total_n_dofs_at_node );
      CHKERRABORT(system.comm().get(),ierr);
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

      unsigned int total_n_dofs_at_elem = 0;

      // We are assuming variables are also numbered 0 to n_vars()-1
      for( unsigned int v = 0; v < system.n_vars(); v++ )
        {
          unsigned int n_dofs_at_elem = elem.n_dofs(system.number(), v);

          if( n_dofs_at_elem > 0 )
            {
              ierr = PetscSectionSetFieldDof( section, local_elem_id, v, n_dofs_at_elem );
              CHKERRABORT(system.comm().get(),ierr);

              total_n_dofs_at_elem += n_dofs_at_elem;
            }
        }

      libmesh_assert_equal_to(total_n_dofs_at_elem, elem.n_dofs(system.number()));

      ierr = PetscSectionSetDof( section, local_elem_id, total_n_dofs_at_elem );
      CHKERRABORT(system.comm().get(),ierr);
    }

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
