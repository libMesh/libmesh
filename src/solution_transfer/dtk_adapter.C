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



#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK

// libMesh includes
#include "libmesh/dtk_adapter.h"
#include "libmesh/dtk_evaluator.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"

// Trilinos includes
#include "libmesh/ignore_warnings.h"
#include <DTK_MeshTypes.hpp>
#include <Teuchos_Comm.hpp>
#include "libmesh/restore_warnings.h"

// C++ includes
#include <vector>
#include <numeric>

namespace libMesh
{

DTKAdapter::DTKAdapter(Teuchos::RCP<const Teuchos::Comm<int>> in_comm, EquationSystems & in_es):
  comm(in_comm),
  es(in_es),
  mesh(in_es.get_mesh()),
  dim(mesh.mesh_dimension())
{
  std::set<unsigned int> semi_local_nodes;
  get_semi_local_nodes(semi_local_nodes);

  num_local_nodes = semi_local_nodes.size();

  vertices.resize(num_local_nodes);
  Teuchos::ArrayRCP<double> coordinates(num_local_nodes * dim);

  // Fill in the vertices and coordinates
  {
    unsigned int i = 0;

    for (const auto & id : semi_local_nodes)
      {
        const Node & node = mesh.node_ref(id);

        vertices[i] = node.id();

        for (unsigned int j=0; j<dim; j++)
          coordinates[(j*num_local_nodes) + i] = node(j);

        i++;
      }
  }

  // Currently assuming all elements are the same!
  DataTransferKit::DTK_ElementTopology element_topology = get_element_topology(mesh.elem_ptr(0));
  unsigned int n_nodes_per_elem = mesh.elem_ptr(0)->n_nodes();

  unsigned int n_local_elem = mesh.n_local_elem();

  Teuchos::ArrayRCP<int> elements(n_local_elem);
  Teuchos::ArrayRCP<int> connectivity(n_nodes_per_elem*n_local_elem);

  // Fill in the elements and connectivity
  {
    unsigned int i = 0;

    for (const auto & elem : as_range(mesh.local_elements_begin(),
                                      mesh.local_elements_end()))
      {
        elements[i] = elem->id();

        for (unsigned int j=0; j<n_nodes_per_elem; j++)
          connectivity[(j*n_local_elem)+i] = elem->node_id(j);

        i++;
      }
  }

  Teuchos::ArrayRCP<int> permutation_list(n_nodes_per_elem);
  std::iota(permutation_list.begin(), permutation_list.end(), 0);

  /*
    if (this->processor_id() == 1)
    sleep(1);

    libMesh::out<<"n_nodes_per_elem: "<<n_nodes_per_elem<<std::endl;

    libMesh::out<<"Dim: "<<dim<<std::endl;

    libMesh::err<<"Vertices size: "<<vertices.size()<<std::endl;
    {
    libMesh::err<<this->processor_id()<<" Vertices: ";

    for (std::size_t i=0; i<vertices.size(); i++)
    libMesh::err<<vertices[i]<<" ";

    libMesh::err<<std::endl;
    }

    libMesh::err<<"Coordinates size: "<<coordinates.size()<<std::endl;
    {
    libMesh::err<<this->processor_id()<<" Coordinates: ";

    for (std::size_t i=0; i<coordinates.size(); i++)
    libMesh::err<<coordinates[i]<<" ";

    libMesh::err<<std::endl;
    }

    libMesh::err<<"Connectivity size: "<<connectivity.size()<<std::endl;
    {
    libMesh::err<<this->processor_id()<<" Connectivity: ";

    for (std::size_t i=0; i<connectivity.size(); i++)
    libMesh::err<<connectivity[i]<<" ";

    libMesh::err<<std::endl;
    }

    libMesh::err<<"Permutation_List size: "<<permutation_list.size()<<std::endl;
    {
    libMesh::err<<this->processor_id()<<" Permutation_List: ";

    for (std::size_t i=0; i<permutation_list.size(); i++)
    libMesh::err<<permutation_list[i]<<" ";

    libMesh::err<<std::endl;
    }

  */
  Teuchos::RCP<MeshContainerType>
    mesh_container = Teuchos::rcp(new MeshContainerType(dim,
                                                        vertices,
                                                        coordinates,
                                                        element_topology,
                                                        n_nodes_per_elem,
                                                        elements,
                                                        connectivity,
                                                        permutation_list));

  // We only have 1 element topology in this grid so we make just one mesh block
  Teuchos::ArrayRCP<Teuchos::RCP<MeshContainerType>> mesh_blocks(1);
  mesh_blocks[0] = mesh_container;

  // Create the MeshManager
  mesh_manager = Teuchos::rcp(new DataTransferKit::MeshManager<MeshContainerType>(mesh_blocks, comm, dim) );

  // Pack the coordinates into a field, this will be the positions we'll ask for other systems fields at
  target_coords = Teuchos::rcp(new DataTransferKit::FieldManager<MeshContainerType>(mesh_container, comm));
}

DTKAdapter::RCP_Evaluator
DTKAdapter::get_variable_evaluator(std::string var_name)
{
  // We try emplacing a nullptr into the "evaluators" map.
  auto [it, emplaced] = evaluators.emplace(var_name, nullptr);

  // If the emplace succeeded, that means it was a new entry in the
  // map, so we need to actually construct the object.
  if (emplaced)
    {
      System * sys = find_sys(var_name);
      it->second = Teuchos::rcp(new DTKEvaluator(*sys, var_name));
    }

  return it->second;
}

Teuchos::RCP<DataTransferKit::FieldManager<DTKAdapter::FieldContainerType>>
DTKAdapter::get_values_to_fill(std::string var_name)
{
  // We try emplacing a nullptr into the "values_to_fill" map.
  auto [it, emplaced] = values_to_fill.emplace(var_name, nullptr);

  // If the emplace succeeded, that means it was a new entry in the
  // map, so we need to actually construct the object.
  if (emplaced)
    {
      Teuchos::ArrayRCP<double> data_space(num_local_nodes);
      Teuchos::RCP<FieldContainerType> field_container = Teuchos::rcp(new FieldContainerType(data_space, 1));
      it->second = Teuchos::rcp(new DataTransferKit::FieldManager<FieldContainerType>(field_container, comm));
    }

  return it->second;
}

void
DTKAdapter::update_variable_values(std::string var_name)
{
  System * sys = find_sys(var_name);
  unsigned int var_num = sys->variable_number(var_name);

  Teuchos::RCP<FieldContainerType> values = values_to_fill[var_name]->field();

  unsigned int i=0;
  // Loop over the values (one for each node) and assign the value of this variable at each node
  for (const auto & value : *values)
    {
      unsigned int node_num = vertices[i];
      const Node & node = mesh.node_ref(node_num);

      if (node.processor_id() == sys->processor_id())
        {
          // The 0 is for the component... this only works for LAGRANGE!
          dof_id_type dof = node.dof_number(sys->number(), var_num, 0);
          sys->solution->set(dof, value);
        }

      i++;
    }

  sys->solution->close();
}


/**
 * Small helper function for finding the system containing the variable.
 *
 * Note that this implies that variable names are unique across all systems!
 */
System *
DTKAdapter::find_sys(std::string var_name)
{
  System * sys = nullptr;

  // Find the system this variable is from
  for (auto i : make_range(es.n_systems()))
    {
      if (es.get_system(i).has_variable(var_name))
        {
          sys = &es.get_system(i);
          break;
        }
    }

  libmesh_assert(sys);

  return sys;
}

DataTransferKit::DTK_ElementTopology
DTKAdapter::get_element_topology(const Elem * elem)
{
  ElemType type = elem->type();

  if (type == EDGE2)
    return DataTransferKit::DTK_LINE_SEGMENT;
  else if (type == TRI3)
    return DataTransferKit::DTK_TRIANGLE;
  else if (type == QUAD4)
    return DataTransferKit::DTK_QUADRILATERAL;
  else if (type == TET4)
    return DataTransferKit::DTK_TETRAHEDRON;
  else if (type == HEX8)
    return DataTransferKit::DTK_HEXAHEDRON;
  else if (type == PYRAMID5)
    return DataTransferKit::DTK_PYRAMID;

  libmesh_error_msg("Element type not supported by DTK!");
}

void
DTKAdapter::get_semi_local_nodes(std::set<unsigned int> & semi_local_nodes)
{
  for (const auto & elem : as_range(mesh.local_elements_begin(), mesh.local_elements_end()))
    for (const Node & node : elem->node_ref_range())
      semi_local_nodes.insert(node.id());
}

} // namespace libMesh

#endif // #ifdef LIBMESH_TRILINOS_HAVE_DTK
