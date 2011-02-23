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

#include "rb_ref_comp_eval.h"
#include "rb_reference_component.h"
#include "libmesh_logging.h"

namespace libMesh
{

RBReferenceComponentEvaluation::RBReferenceComponentEvaluation(RBSystem& rb_sys)
  :
  Parent(rb_sys)
{}

void RBReferenceComponentEvaluation::initialize()
{  
  Parent::initialize();

  RBReferenceComponent& ref_sys = libmesh_cast_ref<RBReferenceComponent&>(rb_sys);
  const unsigned int Nmax   = ref_sys.get_Nmax();
  const unsigned int n_dofs = ref_sys.n_local_interface_functions();
  const unsigned int Q_a    = ref_sys.get_Q_a();

  Aq_g_g.resize(Q_a);
  for(unsigned int i=0; i<Q_a; i++)
  {
    Aq_g_g[i].resize(n_dofs);
    for(unsigned int j=0; j<n_dofs; j++)
    {
      Aq_g_g[i][j].resize(n_dofs);
    }
  }

  Aq_bubble_g.resize(Q_a);
  for(unsigned int i=0; i<Q_a; i++)
  {
    Aq_bubble_g[i].resize(n_dofs);
    for(unsigned int j=0; j<n_dofs; j++)
    {
      Aq_bubble_g[i][j].resize(Nmax);
    }
  }
}

void RBReferenceComponentEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBReferenceComponentEvaluation");

  Parent::write_offline_data_to_files(directory_name);
  
  RBReferenceComponent& ref_sys = libmesh_cast_ref<RBReferenceComponent&>(rb_sys);
  const unsigned int n_dofs = ref_sys.n_local_interface_functions();
  const unsigned int Q_a    = ref_sys.get_Q_a();
  const unsigned int n_bfs  = get_n_basis_functions();

  libmesh_assert( n_bfs <= ref_sys.get_Nmax() );

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // write out Aq_g_g data
    std::ofstream Aq_g_g_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_g_g.dat";
      Aq_g_g_out.open(file_name.str().c_str());
    }
    if ( !Aq_g_g_out.good() )
    {
      libMesh::err << "Error opening Aq_g_g.dat" << std::endl;
      libmesh_error();
    }
    Aq_g_g_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      for(unsigned int i=0; i<n_dofs; i++)
      {
        for(unsigned int j=0; j<n_dofs; j++)
        {
          Aq_g_g_out << std::scientific << Aq_g_g[q_a][i][j] << " ";
        }
      }
    }
    Aq_g_g_out.close();

    // write out Aq_bubble_g data
    std::ofstream Aq_bubble_g_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_bubble_g.dat";
      Aq_bubble_g_out.open(file_name.str().c_str());
    }
    if ( !Aq_bubble_g_out.good() )
    {
      libMesh::err << "Error opening Aq_bubble_g.dat" << std::endl;
      libmesh_error();
    }
    Aq_bubble_g_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      for(unsigned int i=0; i<n_dofs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          Aq_bubble_g_out << std::scientific << Aq_bubble_g[q_a][i][j] << " ";
        }
      }
    }
    Aq_bubble_g_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "RBReferenceComponentEvaluation");
}

void RBReferenceComponentEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBReferenceComponentEvaluation");

  Parent::read_offline_data_from_files(directory_name);

  RBReferenceComponent& ref_sys = libmesh_cast_ref<RBReferenceComponent&>(rb_sys);
  const unsigned int n_dofs = ref_sys.n_local_interface_functions();
  const unsigned int Q_a    = ref_sys.get_Q_a();
  const unsigned int n_bfs  = get_n_basis_functions();

  // Read in Aq_g_g
  std::ifstream Aq_g_g_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_g_g.dat";
    Aq_g_g_in.open(file_name.str().c_str());
  }
  if ( !Aq_g_g_in.good() )
  {
    libMesh::err << "Error opening Aq_g_g.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    for(unsigned int n=0; n<n_dofs; n++)
    {
      for(unsigned int i=0; i<n_dofs; i++)
      {
        Aq_g_g_in >> Aq_g_g[q_a][n][i];
      }
    }
  }
  Aq_g_g_in.close();

  // Read in Aq_bubble_g
  std::ifstream Aq_bubble_g_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_bubble_g.dat";
    Aq_bubble_g_in.open(file_name.str().c_str());
  }
  if ( !Aq_bubble_g_in.good() )
  {
    libMesh::err << "Error opening Aq_bubble_g.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    for(unsigned int n=0; n<n_dofs; n++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        Aq_bubble_g_in >> Aq_bubble_g[q_a][n][i];
      }
    }
  }
  Aq_bubble_g_in.close();

  STOP_LOG("read_offline_data_from_files()", "RBReferenceComponentEvaluation");
}

}
