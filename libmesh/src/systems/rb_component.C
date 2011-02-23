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
#include "rb_component.h"
#include "rb_reference_component.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

namespace libMesh
{

RBComponent::RBComponent (RBReferenceComponent& ref_component_in, RBSystem* source_system_in)
  :
  ref_component(ref_component_in),
  global_id(0)
{
  interface_function_global_indices.resize(ref_component.n_local_interface_functions(), 0);

  source_Fq_g.clear();
  source_Aq_bubble_g.clear();

  neighbor_components.resize(ref_component.get_n_ports());
  // intialize neighbors to NULL
  for(unsigned int i=0; i<neighbor_components.size(); i++)
  {
    neighbor_components[i] = NULL;
  }

  set_source_system(source_system_in);
}

RBComponent::~RBComponent ()
{}

void RBComponent::set_global_id(unsigned int global_id_in)
{
  global_id = global_id_in;
}

unsigned int RBComponent::get_global_id() const
{
  return global_id;
}

void RBComponent::set_source_system(RBSystem* source_system_in)
{
  source_system = source_system_in;

  if(source_system)
  {
    source_Fq_g.resize(source_system->get_Q_f());
    for(unsigned int i=0; i<source_system->get_Q_f(); i++)
    {
      source_Fq_g[i].resize(ref_component.n_local_interface_functions());
    }
    
    unsigned int Q_a  = ref_component.get_Q_a();
    unsigned int Nmax = ref_component.get_Nmax();
    source_Aq_bubble_g.resize(Q_a);
    for(unsigned int i=0; i<Q_a; i++)
    {
      source_Aq_bubble_g[i].resize(ref_component.n_local_interface_functions());
      for(unsigned int j=0; j<ref_component.n_local_interface_functions(); j++)
      {
        source_Aq_bubble_g[i][j].resize(Nmax);
      }
    }
  }
}

void RBComponent::set_component_parameters(const std::vector<Real>& params)
{
  libmesh_assert(params.size() == ref_component.get_n_params());

  component_parameters = params;
}

void RBComponent::set_neighbor_component(unsigned int port, RBComponent& neighbor)
{
  // First make sure we have identified a valid port
  libmesh_assert( port < ref_component.get_n_ports() );
  
  // Then set the neighbor
  neighbor_components[port] = &neighbor;
}

void RBComponent::set_neighbor_components(std::vector< RBComponent* >& neighbors)
{
  libmesh_assert(ref_component.get_n_ports() == neighbors.size());
  
  for(unsigned int i=0; i<neighbors.size(); i++)
  {
    set_neighbor_component(i, *neighbors[i]);
  }
}

void RBComponent::assemble_local_truth_matrix_and_rhs(DenseMatrix<Number>& local_matrix,
                                                      DenseVector<Number>& local_rhs)
{
  unsigned int n_dofs = ref_component.n_local_interface_functions();
  local_rhs.resize(n_dofs);

  // Need to set the ref_component's parameters appropriately
  ref_component.set_current_parameters( get_component_parameters() );

  // First, if necessary, compute the truth bubble function associated with the source
  if(source_system)
  {
    source_system->set_current_parameters( get_component_parameters() );
    source_system->truth_solve(-1);

    AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
    temp->init (source_system->n_dofs(), source_system->n_local_dofs(), false, libMeshEnums::PARALLEL);

    for(unsigned int i=0; i<n_dofs; i++)
    {
      for(unsigned int q_f_source=0; q_f_source<source_system->get_Q_f(); q_f_source++)
      {
        // Need to use non-Dirichlet vectors for the inner products here
        local_rhs(i) += source_system->eval_theta_q_f(q_f_source) *
          source_system->get_non_dirichlet_F_q(q_f_source)->dot(ref_component.get_interface_function_representor(i));
      }

      for(unsigned int q_a=0; q_a<ref_component.get_Q_a(); q_a++)
      {
        temp->zero();
        if(!ref_component.low_memory_mode)
        {
          ref_component.get_non_dirichlet_A_q(q_a)->vector_mult(*temp, *source_system->solution);
        }
        else
        {
          ref_component.assemble_Aq_matrix(q_a, ref_component.matrix, false /* apply_dirichlet_bc */);
          ref_component.matrix->vector_mult(*temp, *source_system->solution);
        }

        local_rhs(i) -= ref_component.eval_theta_q_a(q_a) *
          temp->dot(ref_component.get_interface_function_representor(i));
      }
    }

  }

  ref_component.assemble_local_truth_matrix(local_matrix);
}

void RBComponent::assemble_local_RB_matrix_and_rhs(unsigned int N_rb,
                                                   DenseMatrix<Number>& local_matrix,
                                                   DenseVector<Number>& local_rhs)
{
  unsigned int n_dofs = ref_component.n_local_interface_functions();
  local_rhs.resize(n_dofs);

  // Need to set the ref_component's parameters appropriately
  ref_component.set_current_parameters( get_component_parameters() );

  if(source_system)
  {
    source_system->set_current_parameters( get_component_parameters() );
    unsigned int N = std::min(N_rb, source_system->get_n_basis_functions());
    source_system->rb_eval->RB_solve(N);

    if(source_system->store_basis_functions)
    {
      source_system->load_RB_solution();
    }

    for(unsigned int i=0; i<n_dofs; i++)
    {
      for(unsigned int q_f_source=0; q_f_source<source_system->get_Q_f(); q_f_source++)
      {
        local_rhs(i) += source_system->eval_theta_q_f(q_f_source) * source_Fq_g[q_f_source][i];
      }

      for(unsigned int q_a=0; q_a<ref_component.get_Q_a(); q_a++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          local_rhs(i) -= ref_component.eval_theta_q_a(q_a) *
            source_system->rb_eval->RB_solution(j) * source_Aq_bubble_g[q_a][i][j];
        }
      }
    }
  }

  ref_component.assemble_local_RB_matrix(N_rb, local_matrix);
}

void RBComponent::update_extra_source_data()
{
  START_LOG("update_extra_source_data()", "RBComponent");
  
  if(!source_system)
    return;
    
  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  temp->init (source_system->n_dofs(), source_system->n_local_dofs(), false, libMeshEnums::PARALLEL);

  unsigned int n_dofs = ref_component.n_local_interface_functions();
  for(unsigned int i=0; i<n_dofs; i++)
  {
    // Compute the Fq/interface-function cross terms
    for(unsigned int q_f_source=0; q_f_source<source_system->get_Q_f(); q_f_source++)
    {
      // Need to use non-Dirichlet vectors for the inner products here
      source_Fq_g[q_f_source][i] =
        source_system->get_non_dirichlet_F_q(q_f_source)->dot(ref_component.get_interface_function_representor(i));
    }

    // Compute the Aq/interface-function cross terms
    for(unsigned int q_a=0; q_a<ref_component.get_Q_a(); q_a++)
    {
      for(unsigned int j=0; j<source_system->get_n_basis_functions(); j++)
      {
        temp->zero();
        if(!ref_component.low_memory_mode)
        {
          ref_component.get_non_dirichlet_A_q(q_a)->vector_mult(*temp, source_system->get_basis_function(j));
        }
        else
        {
          ref_component.assemble_Aq_matrix(q_a, ref_component.matrix, false /* apply_dirichlet_bc */);
          ref_component.matrix->vector_mult(*temp, source_system->get_basis_function(j));
        }

        source_Aq_bubble_g[q_a][i][j] =
          temp->dot(ref_component.get_interface_function_representor(i));
      }
    }
  }

  STOP_LOG("update_extra_source_data()", "RBComponent");
}

void RBComponent::load_component_solution(DenseVector<Number>& global_solution)
{
  START_LOG("load_local_solution()", "RBComponent");
  
  ref_component.solution->zero();
  
  if(source_system)
  {
    // We assume that source_system->solution has been
    // set to the appopriate truth source bubble function
    ref_component.solution->add(1., *source_system->solution);
  }

  for(unsigned int index=0; index<ref_component.n_local_interface_functions(); index++)
  {
    // We assume that bubble_functions have been
    // set to the appopriate truth source bubble function
    Number scalar = global_solution(interface_function_global_indices[index]);
    ref_component.solution->add(scalar, *ref_component.bubble_functions[index]);
    ref_component.solution->add(scalar, *ref_component.interface_function_representors[index]);
  }
  ref_component.update();

  STOP_LOG("load_local_solution()", "RBComponent");
}

void RBComponent::write_extra_source_data(const std::string& directory_name)
{
  START_LOG("write_extra_source_data()", "RBComponent");
  
  if(!source_system)
    return;
  
  if(libMesh::processor_id() == 0)
  {
    OStringStream full_directory_name;
    full_directory_name << directory_name << "_"
                        << source_system->name()
                        << "_" << ref_component.name();

    // Make a directory to store all the data files
    if( mkdir(full_directory_name.str().c_str(), 0777) == -1)
    {
      libMesh::out << "In RBComponent::write_extra_source_data, directory "
                   << full_directory_name.str() << " already exists, overwriting contents." << std::endl;
    }

    const unsigned int precision_level = 14;

    // write out source_Fq_g data
    std::ofstream source_Fq_g_out;
    {
      OStringStream file_name;
      file_name << full_directory_name.str() << "/source_Fq_g.dat";
      source_Fq_g_out.open(file_name.str().c_str());
    }
    if ( !source_Fq_g_out.good() )
    {
      libMesh::err << "Error opening source_Fq_g.dat" << std::endl;
      libmesh_error();
    }
    source_Fq_g_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<source_system->get_Q_f(); q_f++)
    {
      for(unsigned int i=0; i<ref_component.n_local_interface_functions(); i++)
      {
        source_Fq_g_out << std::scientific << source_Fq_g[q_f][i] << " ";
      }
    }
    source_Fq_g_out.close();
    
    // write out source_Aq_bubble_g data
    std::ofstream source_Aq_bubble_g_out;
    {
      OStringStream file_name;
      file_name << full_directory_name.str() << "/source_Aq_bubble_g.dat";
      source_Aq_bubble_g_out.open(file_name.str().c_str());
    }
    if ( !source_Aq_bubble_g_out.good() )
    {
      libMesh::err << "Error opening source_Aq_bubble_g.dat" << std::endl;
      libmesh_error();
    }
    source_Aq_bubble_g_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<ref_component.get_Q_a(); q_a++)
    {
      for(unsigned int n=0; n<ref_component.n_local_interface_functions(); n++)
      {
        for(unsigned int i=0; i<source_system->get_n_basis_functions(); i++)
        {
          source_Aq_bubble_g_out << std::scientific << source_Aq_bubble_g[q_a][n][i] << " ";
        }
      }
    }
    source_Aq_bubble_g_out.close();
  }

  STOP_LOG("write_extra_source_data()", "RBComponent");
}

void RBComponent::read_extra_source_data(const std::string& directory_name)
{
  START_LOG("read_extra_source_data()", "RBComponent");

  if(!source_system)
    return;

  OStringStream full_directory_name;
  full_directory_name << directory_name << "_" << source_system->name()
                      << "_" << ref_component.name();

  // Read in source_Fq_g
  std::ifstream source_Fq_g_in;
  {
    OStringStream file_name;
    file_name << full_directory_name.str() << "/source_Fq_g.dat";
    source_Fq_g_in.open(file_name.str().c_str());
  }
  if ( !source_Fq_g_in.good() )
  {
    libMesh::err << "Error opening source_Fq_g.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<source_system->get_Q_f(); q_f++)
  {
    for(unsigned int i=0; i<ref_component.n_local_interface_functions(); i++)
    {
      source_Fq_g_in >> source_Fq_g[q_f][i];
    }
  }
  source_Fq_g_in.close();
  
  // Read in source_Aq_bubble_g
  std::ifstream source_Aq_bubble_g_in;
  {
    OStringStream file_name;
    file_name << full_directory_name.str() << "/source_Aq_bubble_g.dat";
    source_Aq_bubble_g_in.open(file_name.str().c_str());
  }
  if ( !source_Aq_bubble_g_in.good() )
  {
    libMesh::err << "Error opening source_Aq_bubble_g.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<ref_component.get_Q_a(); q_a++)
  {
    for(unsigned int n=0; n<ref_component.n_local_interface_functions(); n++)
    {
      for(unsigned int i=0; i<source_system->get_n_basis_functions(); i++)
      {
        source_Aq_bubble_g_in >> source_Aq_bubble_g[q_a][n][i];
      }
    }
  }
  source_Aq_bubble_g_in.close();

  STOP_LOG("read_extra_source_data()", "RBComponent");
}

}