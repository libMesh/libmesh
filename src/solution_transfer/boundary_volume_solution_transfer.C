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

#include "libmesh/boundary_volume_solution_transfer.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"

namespace libMesh {

void BoundaryVolumeSolutionTransfer::transfer(const Variable & from_var,
                                              const Variable & to_var)
{
  // Determine which direction the transfer is in
  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  unsigned int
    from_dimension = from_sys->get_mesh().mesh_dimension(),
    to_dimension = to_sys->get_mesh().mesh_dimension();

  // Sanity check
  if (from_dimension == to_dimension)
    libmesh_error_msg("Error: Transfer must be from volume mesh to its boundary or vice-versa!");

  if (from_dimension > to_dimension)
    this->transfer_volume_boundary(from_var, to_var);
  else
    this->transfer_boundary_volume(from_var, to_var);
}



void BoundaryVolumeSolutionTransfer::
transfer_volume_boundary(const Variable & from_var, const Variable & to_var)
{
  // Get references to the Systems from the Variables
  System * from_sys = from_var.system(); // volume system
  System * to_sys = to_var.system();     // boundary system

  // Get reference to the BoundaryMesh.  Note: we always loop over the
  // BoundaryMesh since, by definition, it has fewer dofs than the
  // volume mesh.
  const MeshBase & to_mesh = to_sys->get_mesh();

  // Get system number and variable numbers
  const unsigned short int from_sys_number = from_sys->number();
  const unsigned short int from_var_number = from_var.number();

  const unsigned short int to_sys_number = to_sys->number();
  const unsigned short int to_var_number = to_var.number();

  // Get a constant reference to variables, get their number of components
  const unsigned short int from_n_comp = from_var.n_components();
  const unsigned short int to_n_comp = to_var.n_components();

  // Sanity check that the variables have the same number of components
  libmesh_assert_equal_to(from_n_comp, to_n_comp);

  // Construct map from "from" dofs to "to" dofs.
  typedef std::map<numeric_index_type, numeric_index_type> DofMapping;
  DofMapping dof_mapping;

  // Loop through all boundary elements.
  for (const auto & to_elem : to_mesh.active_local_element_ptr_range())
    {
      const Elem * from_elem = to_elem->interior_parent();

      if (!from_elem)
        libmesh_error_msg("Error, transfer must be between a Mesh and its associated BoundaryMesh.");

      // loop through all nodes in each boundary element.
      for (unsigned int node=0; node < to_elem->n_nodes(); node++)
        {
          // Node in boundary element.
          const Node * to_node = to_elem->node_ptr(node);

          for (unsigned int node_id=0; node_id < from_elem->n_nodes(); node_id++)
            {
              // Nodes in interior_parent element.
              const Node * from_node = from_elem->node_ptr(node_id);

              const dof_id_type from_dof = from_node->dof_number(from_sys_number,
                                                                 from_var_number,
                                                                 from_n_comp - 1);

              // See if we've already encountered this DOF in the loop
              // over boundary elements.
              DofMapping::iterator it = dof_mapping.find(from_dof);

              // If we've already mapped this dof, we don't need to map
              // it again or do floating point comparisons.
              if (it == dof_mapping.end() &&
                  from_node->absolute_fuzzy_equals(*to_node, TOLERANCE))
                {
                  // Global dof_index for node in BoundaryMesh.
                  const dof_id_type to_dof = to_node->dof_number(to_sys_number,
                                                                 to_var_number,
                                                                 to_n_comp - 1);

                  // Keep track of the volume system dof index which is needed.
                  dof_mapping[from_dof] = to_dof;
                }
            }
        }
    }

  // Construct a vector of the indices needed from the Volume system's
  // global solution vector on this processor.
  std::vector<numeric_index_type> needed_indices;
  needed_indices.reserve(dof_mapping.size());
  {
    DofMapping::iterator
      it = dof_mapping.begin(),
      end = dof_mapping.end();

    for (; it!=end; ++it)
      needed_indices.push_back(it->first);
  }

  // Communicate just the required values without making a copy of the
  // global solution vector.
  std::vector<Number> needed_values;
  from_sys->solution->localize(needed_values, needed_indices);

  // Loop over DofMapping again, assigning values in the
  // Boundary System solution vector.
  {
    DofMapping::iterator
      it = dof_mapping.begin(),
      end = dof_mapping.end();

    for (unsigned idx=0; it!=end; ++it, ++idx)
      to_sys->solution->set(it->second, needed_values[idx]);
  }
}


void BoundaryVolumeSolutionTransfer::
transfer_boundary_volume(const Variable & from_var, const Variable & to_var)
{
  // Get references to the Systems from the Variables
  System * from_sys = from_var.system(); // boundary system
  System * to_sys = to_var.system();     // volume system

  // Get reference to BoundaryMesh.
  const MeshBase & from_mesh = from_sys->get_mesh();

  // DofMap for BoundaryMesh
  const DofMap & dof_map = from_sys->get_dof_map();

  // Get system number and variable numbers
  const unsigned short int to_sys_number = to_sys->number();
  const unsigned short int to_var_number = to_var.number();
  const unsigned short int from_var_number = from_var.number();
  const unsigned short int to_n_comp = to_var.n_components();

  // In order to get solution vectors from BoundaryMesh
  std::vector<dof_id_type> from_dof_indices;
  std::vector<Number> value;

  // Loop through all boundary elements.
  for (const auto & from_elem : from_mesh.active_local_element_ptr_range())
    {
      const Elem * to_elem = from_elem->interior_parent();

      if (!to_elem)
        libmesh_error_msg("Error, transfer must be between a Mesh and its associated BoundaryMesh.");

      // Get dof indices for phi2 for all nodes on this boundary element
      dof_map.dof_indices(from_elem, from_dof_indices, from_var_number);

      // Get values from BoundaryMesh for this element
      from_sys->current_local_solution->get(from_dof_indices, value);

      // loop through all nodes in each boundary element.
      for (unsigned int node=0; node < from_elem->n_nodes(); node++)
        {
          // Node in boundary element.
          const Node * from_node = from_elem->node_ptr(node);

          for (unsigned int node_id=0; node_id < to_elem->n_nodes(); node_id++)
            {
              // Nodes in interior_parent element.
              const Node * to_node = to_elem->node_ptr(node_id);

              // Match BoundaryNode & VolumeNode.
              if (to_node->absolute_fuzzy_equals(*from_node, TOLERANCE))
                {
                  // Global dof_index for node in VolumeMesh.
                  const dof_id_type to_dof = to_node->dof_number(to_sys_number,
                                                                 to_var_number,
                                                                 to_n_comp - 1);

                  // Assign values to boundary in VolumeMesh.
                  to_sys->solution->set(to_dof, value[node]);
                }
            }
        }
    }
}

} // namespace libMesh
