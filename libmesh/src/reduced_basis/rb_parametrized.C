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

#include "rb_parametrized.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBParameters implementation

RBParametrized::RBParametrized()
  :
  parameters_initialized(false)
{
  libmesh_experimental();
  
  current_parameters.clear();
  mu_min_vector.clear();
  mu_max_vector.clear();
}

RBParametrized::~RBParametrized()
{
  this->clear();
}

void RBParametrized::clear()
{
  current_parameters.clear();
  mu_min_vector.clear();
  mu_max_vector.clear();
  parameters_initialized = false;
}

void RBParametrized::initialize_parameters(std::vector<Real> mu_min_in,
                                           std::vector<Real> mu_max_in,
                                           std::vector<Real> mu_in)
{
  // Check that the min/max vectors are valid
  {
    const std::string err_string = "Error: Invalid mu_min/mu_max in RBParameters constructor.";
    bool valid_min_max = (mu_min_in.size() == mu_max_in.size());
    if(!valid_min_max)
    {
      libMesh::err << err_string << std::endl;
    }
    else
    {
      for(unsigned int i=0; i<mu_min_in.size(); i++)
      {
        if(mu_min_in[i] > mu_max_in[i])
        {
          libMesh::err << err_string << std::endl;
        }
      }
    }
  }
  
  mu_min_vector = mu_min_in;
  mu_max_vector = mu_max_in;
  // Need to resize so that error checking in set_current_parameters works properly
  current_parameters.resize(mu_min_vector.size());

  parameters_initialized = true;
  
  set_current_parameters(mu_in);
}

void RBParametrized::initialize_parameters(RBParametrized& rb_parametrized)
{
  const unsigned int n_params = rb_parametrized.get_n_params();
  std::vector<Real> parameters_min_vector(n_params);
  std::vector<Real> parameters_max_vector(n_params);
  for(unsigned int i=0; i<n_params; i++)
  {
    parameters_min_vector[i] = rb_parametrized.get_parameter_min(i);
    parameters_max_vector[i] = rb_parametrized.get_parameter_max(i);
  }

  initialize_parameters(parameters_min_vector,
                        parameters_max_vector,
                        rb_parametrized.get_current_parameters());
}

unsigned int RBParametrized::get_n_params() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_n_params" << std::endl;
    libmesh_error();
  }

  return current_parameters.size();
}

void RBParametrized::set_current_parameters(const std::vector<Real>& params)
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::set_current_parameters" << std::endl;
    libmesh_error();
  }

  if(!valid_params(params))
  {
    libMesh::err << "Invalid input parameters in RBParametrized::set_current_parameters" << std::endl;

    for(unsigned int i=0; i<params.size(); i++)
    {
      libMesh::out << "params i = " << params[i] << ", min = " << get_parameter_min(i) << ", max = " << get_parameter_max(i) << std::endl;
    }

    libmesh_error();
  }
  else
  {
    current_parameters = params;
  }
}

std::vector<Real> RBParametrized::get_current_parameters() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_current_parameters" << std::endl;
    libmesh_error();
  }

  return current_parameters;
}

Real RBParametrized::get_parameter_min(unsigned int i) const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameter_min" << std::endl;
    libmesh_error();
  }

  if(i >= get_n_params())
  {
    libMesh::err << "Error: index too large in RBParametrized::get_parameter_min" << std::endl;
    libmesh_error();
  }

  return mu_min_vector[i];
}

Real RBParametrized::get_parameter_max(unsigned int i) const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameter_max" << std::endl;
    libmesh_error();
  }

  if(i >= get_n_params())
  {
    libMesh::err << "Error: index too large in RBParametrized::get_parameter_min" << std::endl;
    libmesh_error();
  }

  return mu_max_vector[i];
}

void RBParametrized::print_current_parameters()
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::print_current_parameters" << std::endl;
    libmesh_error();
  }

  for(unsigned int j=0; j<get_n_params(); j++)
  {
    libMesh::out << "mu[" << j << "] = " << current_parameters[j] << std::endl;
  }
  libMesh::out << std::endl;
}

bool RBParametrized::valid_params(const std::vector<Real>& params)
{
  bool valid = ( params.size() == get_n_params() );

  if(!valid)
  {
   return false;
  }
  else
  {
    for(unsigned int i=0; i<params.size(); i++)
    {
      valid = valid && ( (get_parameter_min(i) <= params[i]) &&
                         (params[i] <= get_parameter_max(i)) );
    }
  }

  return valid;
}

} // namespace libMesh
