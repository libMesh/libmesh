// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "rb_component_system.h"
#include "rb_reference_component.h"
#include "libmesh_logging.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

namespace libMesh
{

RBComponentSystem::RBComponentSystem ()
  :
  n_global_dofs(0)
{}

RBComponentSystem::~RBComponentSystem ()
{
  const_component_iterator iter = components.begin();
  for( ; iter != components.end(); iter++)
  {
    RBComponent* comp = *iter;
    if(comp)
    {
      delete comp;
      comp = NULL;
    }
  }
  components.clear();
}

unsigned int RBComponentSystem::get_n_global_dofs() const
{
  return n_global_dofs;
}

RBComponent& RBComponentSystem::add_component(RBReferenceComponent& ref_component, RBSystem* source_system)
{
  RBComponent* new_comp = new RBComponent(ref_component, source_system);
  new_comp->set_global_id(components.size());
  
  components.push_back( new_comp );
  
  return *new_comp;
}

void RBComponentSystem::initialize_component_system(bool do_not_assemble)
{
  // Initialize the reference components
  for(unsigned int i=0; i<reference_components.size(); i++)
  {
    reference_components[i]->initialize_RB_system(do_not_assemble);
  }
  
  // Initialize the source systems
  for(unsigned int i=0; i<source_systems.size(); i++)
  {
    source_systems[i]->initialize_RB_system(do_not_assemble);
  }
}

void RBComponentSystem::train_and_write_all_components(const std::string& directory_name)
{
  // Train the reference components
  for(unsigned int i=0; i<reference_components.size(); i++)
  {
    reference_components[i]->train_and_write_all_reduced_bases(directory_name);
  }

  // Train the source systems
  for(unsigned int i=0; i<source_systems.size(); i++)
  {
    source_systems[i]->train_reduced_basis();

    OStringStream subdirectory_name;
    subdirectory_name << directory_name << "_" << source_systems[i]->name();
    source_systems[i]->write_offline_data_to_files(subdirectory_name.str());
  }
    
  // Write out any extra data from the source systems
  // for each physical component
  for(unsigned int i=0; i<n_components(); i++)
  {
    RBComponent& comp = get_component(i);
      
    if(comp.source_system)
    {
      comp.update_extra_source_data();
      comp.write_extra_source_data(directory_name);
    }
  }

}

void RBComponentSystem::read_in_all_components(const std::string& directory_name)
{
  libMesh::out << "Read in all the offline data for the component system..." << std::endl;
  
  // Read in the reference components
  for(unsigned int i=0; i<reference_components.size(); i++)
  {
    reference_components[i]->read_in_all_reduced_bases(directory_name);
  }

  for(unsigned int i=0; i<source_systems.size(); i++)
  {
    OStringStream source_directory_name;
    source_directory_name << directory_name << "_" << source_systems[i]->name();
    source_systems[i]->read_offline_data_from_files(source_directory_name.str());
  }

  // Read in any extra data from the source systems
  // for each physical component
  for(unsigned int i=0; i<n_components(); i++)
  {
    RBComponent& comp = get_component(i);
      
    if(comp.source_system)
    {
      comp.read_extra_source_data(directory_name);
    }
  }
  
  libMesh::out << "Finished reading offline data for the component system." << std::endl;
}

RBComponent& RBComponentSystem::get_component(unsigned int id)
{
  libmesh_assert( id < n_components() );
  
  return *components[id];
}

void RBComponentSystem::renumber_interface_functions()
{
  unsigned int next_available_index = 0;

  for(unsigned int i=0; i<components.size(); i++)
  {
          
    RBComponent& comp = *components[i];
    
    // We may need to identifity which port isn't "fresh",
    // if so, store the value in unfresh_port_num
    unsigned int unfresh_port_num = -1;
    
    // Loop over the interface functions
    for(unsigned int local_dof=0; local_dof<comp.ref_component.n_local_interface_functions(); local_dof++)
    {
      // what ports does this interface function have support on?
      std::vector<unsigned int> ports_with_support = comp.ref_component.get_ports_with_support(local_dof);

      // fresh_port indicates if we have encountered this dof (i.e. interface function)
      // before in looping over the components
      bool fresh_dof = true;

      for(unsigned int index=0; index<ports_with_support.size(); index++)
      {
        unsigned int port_num = ports_with_support[index];
        RBComponent* neighbor = comp.get_neighbor_component(port_num);
      
        // neighbor_exists indicates if we have a neighbor on port j
        bool neighbor_exists = (neighbor!=NULL);
      
        // fresh_port indicates if we have encountered port j before.
        // if there is no neighbor on port j, then it must be a fresh port
        bool fresh_port = !neighbor_exists;
      
        // if there is a neighbor, then check component IDs to find out
        // if we've already looped over the neighbor component
        if(neighbor_exists)
        {
          libmesh_assert(neighbor->get_global_id() != comp.get_global_id());
        
          fresh_port = (neighbor->get_global_id() > comp.get_global_id());
        }
        
        // If any of the ports aren't "fresh", then set fresh_dof to false
        // Otherwise we retain the default value fresh_dof = true.
        if(!fresh_port)
        {
          unfresh_port_num = port_num;
          fresh_dof = false;
          break;
        }
      }
      
      if(fresh_dof)
      {
        comp.interface_function_global_indices[local_dof] = next_available_index;

        next_available_index++;
      }
      else // we have already encountered this local dof from a neighboring element
      {
        RBComponent* neighbor = comp.get_neighbor_component(unfresh_port_num);
        libmesh_assert( neighbor != NULL );

        // find out which of neighbor's ports connects to comp
        for(unsigned int port=0; port<neighbor->ref_component.get_n_ports(); port++)
        {
          if(neighbor->get_neighbor_component(port) != NULL)
          {
            if(neighbor->get_neighbor_component(port)->get_global_id() == comp.get_global_id())
            {
              // We have found neighbor's port that connects to comp.
              // Now copy the global interface function indices from neighbor
              // over to comp.
              // This assumes that the interface functions are (locally) numbered
              // on the common port in the same order on comp and neighbor.
              //
              // This condition can be a bit tricky for components with support
              // on multiple interfaces. For example, we can't use the dof numbering:
              //   3----2
              //   |    |
              //   |    |
              //   0----1
              // where each dof has support on two interfaces ("bilinear").
              // Instead, use:
              //   2----3
              //   |    |
              //   |    |
              //   0----1
              // and don't rotate the elements. This way, dofs are always added
              // in the same sequence on each side of the "mesh".
              //
              // TODO: A more general way to proceed here would be to define
              // a Port class, which has "internal" dofs, but can also store
              // nodes, so that reference dofs can be associated with physical
              // nodes and therefore avoid this "rotation" issue.

              unsigned int n_dofs_on_port = comp.ref_component.get_local_dofs_on_port(unfresh_port_num).size();
              libmesh_assert(n_dofs_on_port == neighbor->ref_component.get_local_dofs_on_port(port).size());

              for(unsigned int k=0; k<n_dofs_on_port; k++)
              {
                unsigned int comp_local_dof = comp.ref_component.get_local_dofs_on_port(unfresh_port_num)[k];
                
                if(comp_local_dof == local_dof)
                {
                  unsigned int neighbor_local_dof = neighbor->ref_component.get_local_dofs_on_port(port)[k];
        
                  comp.interface_function_global_indices[comp_local_dof] =
                    neighbor->interface_function_global_indices[neighbor_local_dof];
                }
              }

              break;
            }
          }
        }

      }
    }

  }
  
  n_global_dofs = next_available_index;
}

void RBComponentSystem::truth_solve()
{
  START_LOG("truth_solve()", "RBComponentSystem");

  assemble_truth_global_matrix_and_rhs();

  // Note that the lu_solve call changes global_matrix
  global_matrix.lu_solve(global_rhs, global_solution);

  STOP_LOG("truth_solve()", "RBComponentSystem");
}

void RBComponentSystem::RB_solve(unsigned int N_rb)
{
  START_LOG("RB_solve()", "RBComponentSystem");
  
  assemble_RB_global_matrix_and_rhs(N_rb);

  // Note that the lu_solve call changes global_matrix
  global_matrix.lu_solve(global_rhs, global_solution);
  
  STOP_LOG("RB_solve()", "RBComponentSystem");
}

void RBComponentSystem::assemble_truth_global_matrix_and_rhs()
{
  // resize the global matrix and vector
  unsigned int size = get_n_global_dofs();
  global_matrix.resize(size, size);
  global_rhs.resize   (size);

  // Also, create a local matrix and vector in which to assemble
  // the contribution from each component
  DenseMatrix<Number> local_matrix;
  DenseVector<Number> local_rhs;

  const_component_iterator iter = components.begin();
  for( ; iter != components.end(); iter++)
  {
    RBComponent* comp = *iter;
    comp->assemble_local_truth_matrix_and_rhs(local_matrix, local_rhs);

    add_local_matrix(local_matrix, local_rhs, comp->interface_function_global_indices);
  }
}

void RBComponentSystem::assemble_RB_global_matrix_and_rhs(unsigned int N_rb)
{
  // resize the global matrix and vector
  unsigned int size = get_n_global_dofs();
  global_matrix.resize(size, size);
  global_rhs.resize   (size);

  // Also, create a local matrix and vector in which to assemble
  // the contribution from each component
  DenseMatrix<Number> local_matrix;
  DenseVector<Number> local_rhs;

  const_component_iterator iter = components.begin();
  for( ; iter != components.end(); iter++)
  {
    RBComponent* comp = *iter;
    comp->assemble_local_RB_matrix_and_rhs(N_rb, local_matrix, local_rhs);
    
    add_local_matrix(local_matrix, local_rhs, comp->interface_function_global_indices);
  }
}

void RBComponentSystem::add_local_matrix(const DenseMatrix<Number>& local_matrix,
                                         const DenseVector<Number>& local_rhs,
                                         std::vector<unsigned int>& global_indices)
{
  unsigned int size = global_indices.size();
  
  libmesh_assert( local_matrix.m() == size );
  libmesh_assert( local_matrix.n() == size );
  libmesh_assert( local_rhs.size() == size );
  
  for(unsigned int i=0; i<size; i++)
  {
    unsigned int global_row = global_indices[i];
    global_rhs(global_row) += local_rhs(i);

    for(unsigned int j=0; j<size; j++)
    {
      unsigned int global_col = global_indices[j];
      global_matrix(global_row,global_col) += local_matrix(i,j);
    }
  }
}

} // namespace libMesh