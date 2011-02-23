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

// rbOOmit includes
#include "rb_reference_component.h"
#include "rb_ref_comp_eval.h"

// LibMesh includes
#include "numeric_vector.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "libmesh_logging.h"

namespace libMesh
{

RBReferenceComponent::RBReferenceComponent(EquationSystems& es,
                         const std::string& name,
                         const unsigned int number)
                         : Parent(es,name,number)
{
}

RBReferenceComponent::~RBReferenceComponent()
{
  this->clear();
}

void RBReferenceComponent::clear()
{
  Parent::clear();
  
  // Clear the interface functions
  for(unsigned int i=0; i < interface_function_representors.size(); i++)
  {
    if( interface_function_representors[i] )
    {
      interface_function_representors[i]->clear();
      delete interface_function_representors[i];
    }
  }
  interface_function_representors.clear();
  
  // Clear local_dofs_on_port
  for(unsigned int i=0; i<local_dofs_on_port.size(); i++)
  {
    local_dofs_on_port[i].clear();
  }
  local_dofs_on_port.clear();

  // Clear ports_with_support
  for(unsigned int i=0; i<ports_with_support.size(); i++)
  {
    ports_with_support[i].clear();
  }
  ports_with_support.clear();

  // Clear the truth bubble functions
  for(unsigned int i=0; i<bubble_functions.size(); i++)
  {
    if(bubble_functions[i])
    {
      bubble_functions[i]->clear();
      delete bubble_functions[i];
    }
  }
  bubble_functions.clear();
}

void RBReferenceComponent::add_new_rb_evaluation_object()
{
  RBReferenceComponentEvaluation* e =
    new RBReferenceComponentEvaluation(*this);
  rb_evaluation_objects.push_back(e);
}

void RBReferenceComponent::set_n_ports(unsigned int n_ports_in)
{
  n_ports = n_ports_in;

  local_dofs_on_port.resize(n_ports);
}

std::vector<unsigned int> RBReferenceComponent::get_local_dofs_on_port(unsigned int i) const
{
  libmesh_assert(i < get_n_ports());

  return local_dofs_on_port[i];
}

std::vector<unsigned int> RBReferenceComponent::get_ports_with_support(unsigned int i) const
{
  libmesh_assert(i < n_local_interface_functions());

  return ports_with_support[i];
}

void RBReferenceComponent::add_interface_function(unsigned int port_num,
                                                  Number g(const Point& p,
                                                           const Parameters& ,
                                                           const std::string& ,
                                                           const std::string& ))
{
  std::vector<unsigned int> port_nums(1); port_nums[0] = port_num;
  add_interface_function(port_nums, g);
}

void RBReferenceComponent::add_interface_function(std::vector<unsigned int>& port_nums,
                                                  Number g(const Point& p,
                                                           const Parameters& ,
                                                           const std::string& ,
                                                           const std::string& ))
{
  START_LOG("add_interface_function()", "RBReferenceComponent");
  
  unsigned int local_index = n_local_interface_functions();
  
  for(unsigned int i=0; i<port_nums.size(); i++)
  {
    unsigned int port_num = port_nums[i];
    local_dofs_on_port[port_num].push_back(local_index);
  }
  ports_with_support.push_back( port_nums );

  interface_function_representors.push_back(NULL); // placeholder for representor vector
  
  interface_function_representors[local_index] = NumericVector<Number>::build().release();
  interface_function_representors[local_index]->init (this->n_dofs(), this->n_local_dofs(),
                                                      this->get_dof_map().get_send_list(), false,
                                                      GHOSTED);

  this->project_vector(g,
                       NULL,
                       this->get_equation_systems().parameters,
                       *interface_function_representors[local_index]);

  STOP_LOG("add_interface_function()", "RBReferenceComponent");
}

void RBReferenceComponent::add_interface_function(unsigned int port_num,
                                                  NumericVector<Number>& vec)
{
  std::vector<unsigned int> port_nums(1); port_nums[0] = port_num;
  add_interface_function(port_nums, vec);
}

void RBReferenceComponent::add_interface_function(std::vector<unsigned int>& port_nums,
                                                  NumericVector<Number>& vec)
{
  START_LOG("add_interface_function()", "RBReferenceComponent");

  unsigned int local_index = n_local_interface_functions();

  for(unsigned int i=0; i<port_nums.size(); i++)
  {
    unsigned int port_num = port_nums[i];
    local_dofs_on_port[port_num].push_back(local_index);
  }
  ports_with_support.push_back( port_nums );
  
  interface_function_representors.push_back(NULL); // placeholder for representor vector

  *interface_function_representors[local_index] = vec;

  STOP_LOG("add_interface_function()", "RBReferenceComponent");
}

void RBReferenceComponent::initialize_RB_system (bool online_mode)
{
  Parent::initialize_RB_system(online_mode);

  if(initialize_calN_dependent_data)
  {
    bubble_functions.resize(n_local_interface_functions());
    for(unsigned int index=0; index<bubble_functions.size(); index++)
    {
      bubble_functions[index] = NumericVector<Number>::build().release();
      bubble_functions[index]->init (this->n_dofs(), this->n_local_dofs(),
                                           this->get_dof_map().get_send_list(), false,
                                           GHOSTED);
      
    }
  }

}

void RBReferenceComponent::attach_component_affine_operator(theta_q_fptr theta_q_a,
                                                            affine_assembly_fptr A_q_intrr_assembly,
                                                            affine_assembly_fptr A_q_bndry_assembly)
{
  Parent::attach_A_q(theta_q_a, A_q_intrr_assembly, A_q_bndry_assembly);

  // the F_q assembly functions should not be called directly, hence set them to NULL
  Parent::attach_F_q(theta_q_a, NULL, NULL);
}

void RBReferenceComponent::assemble_local_truth_matrix(DenseMatrix<Number>& local_matrix)
{
  START_LOG("assemble_local_truth_matrix()", "RBReferenceComponent");
  
  unsigned int n_dofs = n_local_interface_functions();
  local_matrix.resize(n_dofs, n_dofs);

  
  // Compute and store the truth bubble functions associated with each interface function
  for(unsigned int index=0; index<n_local_interface_functions(); index++)
  {
    update_F_q_for_interface_function(index);
    
    truth_solve(-1);
    *bubble_functions[index] = *current_local_solution;
  }

  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  
  for(unsigned int i=0; i<n_dofs; i++)
  {
    for(unsigned int j=0; j<n_dofs; j++)
    {
      temp->zero();
      temp->add(1., *bubble_functions[j]);
      temp->add(1., *interface_function_representors[j]);

      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        inner_product_storage_vector->zero();
        if(!low_memory_mode)
        {
          get_non_dirichlet_A_q(q_a)->vector_mult(*inner_product_storage_vector, *temp);
        }
        else
        {
          assemble_Aq_matrix(q_a, matrix, false /* apply_dirichlet_bc */);
          matrix->vector_mult(*inner_product_storage_vector, *temp);
        }

        local_matrix(i,j) += eval_theta_q_a(q_a) *
          inner_product_storage_vector->dot(*interface_function_representors[i]);
      }
    }
  }

  STOP_LOG("assemble_local_truth_matrix()", "RBReferenceComponent");
}

void RBReferenceComponent::assemble_local_RB_matrix(unsigned int N_rb,
                                                    DenseMatrix<Number>& local_matrix)
{
  START_LOG("assemble_local_RB_matrix()", "RBReferenceComponent");
  
  unsigned int n_dofs = n_local_interface_functions();
  local_matrix.resize(n_dofs, n_dofs);

  for(unsigned int interface_function_index=0; interface_function_index<n_dofs; interface_function_index++)
  {    
    // locate the reduced basis data associated with the current interface function
    rb_eval = rb_evaluation_objects[interface_function_index];
    RBReferenceComponentEvaluation* ref_comp_eval =
      libmesh_cast_ptr<RBReferenceComponentEvaluation*>(rb_eval);

    for(unsigned int j=0; j<n_dofs; j++)
    {
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        local_matrix(interface_function_index,j) += eval_theta_q_a(q_a) *
          ref_comp_eval->Aq_g_g[q_a][interface_function_index][j];
      }
    }

    // Solve the system corresponding to the current interface function
    unsigned int N = std::min(N_rb, get_n_basis_functions());
    ref_comp_eval->RB_solve(N);

    if(store_basis_functions)
    {
      load_RB_solution();
      *bubble_functions[interface_function_index] = *current_local_solution;
    }

    for(unsigned int i=0; i<n_dofs; i++)
    {
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          local_matrix(i,interface_function_index) += eval_theta_q_a(q_a) *
            ref_comp_eval->RB_solution(j) * ref_comp_eval->Aq_bubble_g[q_a][i][j];
        }
      }
    }
  }

  STOP_LOG("assemble_local_RB_matrix()", "RBReferenceComponent");
}

void RBReferenceComponent::train_and_write_all_reduced_bases(const std::string& directory_name)
{
  START_LOG("train_and_write_all_reduced_bases()", "RBReferenceComponent");

  // Now perform a Greedy for each interface function and combine the results
  for(unsigned int index=0; index<n_local_interface_functions(); index++)
  {
    // Clear the data and update the F_q vectors according
    // to the current interface_function
    clear_basis_function_dependent_data();
    update_F_q_for_interface_function(index);
    
    // Train a new set of basis functions
    train_reduced_basis();
    
    // Update the data vectors needed to assemble the local matrix
    update_data_for_local_matrix();
    
    OStringStream subdirectory_name;
    subdirectory_name << directory_name << "_" << name() << "_bubble_" << index;
    write_offline_data_to_files( subdirectory_name.str() );
  }
  
  STOP_LOG("train_and_write_all_reduced_bases()", "RBReferenceComponent");
}

void RBReferenceComponent::read_in_all_reduced_bases(const std::string& directory_name)
{
  START_LOG("read_in_all_reduced_bases()", "RBReferenceComponent");

  // Now perform a Greedy for each interface function and combine the results
  for(unsigned int index=0; index<n_local_interface_functions(); index++)
  {
    if(index==0)
    {
      // Clear any basis functions that might already be
      // in rb_eval
      clear_basis_function_dependent_data();
    }
    else
    {
      add_new_rb_evaluation_object();
      
      libmesh_assert( index == (rb_evaluation_objects.size()-1) );
      
      rb_eval = rb_evaluation_objects[index];
      rb_eval->initialize();
    }
    
    OStringStream subdirectory_name;
    subdirectory_name << directory_name << "_" << name() << "_bubble_" << index;
    read_offline_data_from_files( subdirectory_name.str() );
  }
  
  STOP_LOG("read_in_all_reduced_bases()", "RBReferenceComponent");
}

void RBReferenceComponent::update_F_q_for_interface_function(unsigned int index)
{
  START_LOG("update_F_q_for_interface_function()", "RBReferenceComponent");

  if(index >= n_local_interface_functions())
  {
    libMesh::err << "Error: We must have index < n_local_interface_functions() "
                 << "in update_F_q_for_interface_function."
                 << std::endl;
    libmesh_error();
  }

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    libmesh_assert(!is_F_EIM_function(q_f)); // Haven't implemented this for the EIM yet.

    get_F_q(q_f)->zero();
    if(!low_memory_mode)
    {
      // Perform a matrix multiplication with the non-Dirichlet A
      // operator with the specified interface function
      get_non_dirichlet_A_q(q_f)->vector_mult(*get_F_q(q_f),
                                              *interface_function_representors[index]);
    }
    else
    {
      assemble_Aq_matrix(q_f, matrix, false /* apply_dirichlet_bc */);
      matrix->vector_mult(*get_F_q(q_f),
                          *interface_function_representors[index]);
    }
    // Negate the rhs vectors, as required by the bubble function formulation
    get_F_q(q_f)->scale(-1.);
    
    if(store_non_dirichlet_operators)
    {
      *get_non_dirichlet_F_q(q_f) = *get_F_q(q_f);
    }
    
    // Finally, impose the dirichlet BCs on F_q
    // We need to do this directly because we assembled
    // the F_q using an A_q operator with no Dirichlet BCs
    zero_dirichlet_dofs_on_vector(*get_F_q(q_f));
  }

  STOP_LOG("update_F_q_for_interface_function()", "RBReferenceComponent");
}

void RBReferenceComponent::update_data_for_local_matrix()
{
  START_LOG("update_data_for_local_matrix()", "RBReferenceComponent");

  RBReferenceComponentEvaluation* ref_comp_eval =
    libmesh_cast_ptr<RBReferenceComponentEvaluation*>(rb_eval);

  unsigned int n_dofs = n_local_interface_functions();

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int i=0; i<n_dofs; i++)
    {
      for(unsigned int j=0; j<n_dofs; j++)
      {
        inner_product_storage_vector->zero();
        if(!low_memory_mode)
        {
          get_non_dirichlet_A_q(q_a)->vector_mult(*inner_product_storage_vector, *interface_function_representors[j]);
        }
        else
        {
          assemble_Aq_matrix(q_a, matrix, false /* apply_dirichlet_bc */);
          matrix->vector_mult(*inner_product_storage_vector, *interface_function_representors[j]);
        }

        ref_comp_eval->Aq_g_g[q_a][i][j] =
          inner_product_storage_vector->dot(*interface_function_representors[i]);
      }
    }
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int i=0; i<n_dofs; i++)
    {
      for(unsigned int j=0; j<get_n_basis_functions(); j++)
      {
        inner_product_storage_vector->zero();
        if(!low_memory_mode)
        {
          get_non_dirichlet_A_q(q_a)->vector_mult(*inner_product_storage_vector, get_basis_function(j));
        }
        else
        {
          assemble_Aq_matrix(q_a, matrix, false /* apply_dirichlet_bc */);
          matrix->vector_mult(*inner_product_storage_vector, get_basis_function(j));
        }

        ref_comp_eval->Aq_bubble_g[q_a][i][j] =
          inner_product_storage_vector->dot(*interface_function_representors[i]);
      }
    }
  }

  STOP_LOG("update_data_for_local_matrix()", "RBReferenceComponent");
}

}
