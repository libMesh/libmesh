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

// libMesh includes
#include "libmesh/xdr_cxx.h"

// rbOOmit includes
#include "libmesh/rb_parametrized.h"

// C++ includes
#include <sstream>

namespace libMesh
{

// ------------------------------------------------------------
// RBParameters implementation

RBParametrized::RBParametrized()
  :
  verbose_mode(false),
  parameters_initialized(false)
{
  libmesh_experimental();

  parameters.clear();
  parameters_min.clear();
  parameters_max.clear();
}

RBParametrized::~RBParametrized()
{
  this->clear();
}

void RBParametrized::clear()
{
  parameters.clear();
  parameters_min.clear();
  parameters_max.clear();
  parameters_initialized = false;
}

void RBParametrized::initialize_parameters(const RBParameters& mu_min_in,
                                           const RBParameters& mu_max_in,
                                           const RBParameters& mu_in)
{
  // Check that the min/max vectors are valid
  {
    const std::string err_string = "Error: Invalid mu_min/mu_max in RBParameters constructor.";
    bool valid_min_max = (mu_min_in.n_parameters() == mu_max_in.n_parameters());
    if(!valid_min_max)
    {
      libMesh::err << err_string << std::endl;
    }
    else
    {
      RBParameters::const_iterator it     = mu_min_in.begin();
      RBParameters::const_iterator it_end = mu_min_in.end();
      for( ; it != it_end; ++it)
      {
        std::string param_name = it->first;
        if(mu_min_in.get_value(param_name) > mu_max_in.get_value(param_name))
        {
          libMesh::err << err_string << std::endl;
        }
      }
    }
  }

  parameters_min = mu_min_in;
  parameters_max = mu_max_in;

  parameters_initialized = true;
  set_parameters(mu_in);
}

void RBParametrized::initialize_parameters(const RBParametrized& rb_parametrized)
{
  initialize_parameters(rb_parametrized.get_parameters_min(),
                        rb_parametrized.get_parameters_max(),
                        rb_parametrized.get_parameters());
}

unsigned int RBParametrized::get_n_params() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_n_params" << std::endl;
    libmesh_error();
  }

  libmesh_assert_equal_to ( parameters_min.n_parameters(), parameters_max.n_parameters() );

  return parameters_min.n_parameters();
}

void RBParametrized::set_parameters(const RBParameters& params)
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::set_current_parameters" << std::endl;
    libmesh_error();
  }

  valid_params(params); // Terminates if params has the wrong number of parameters

  // Make a copy of params (default assignment operator just does memberwise copy, which is sufficient here)
  this->parameters = params;
}

const RBParameters& RBParametrized::get_parameters() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_current_parameters" << std::endl;
    libmesh_error();
  }

  return parameters;
}

const RBParameters& RBParametrized::get_parameters_min() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameters_min" << std::endl;
    libmesh_error();
  }

  return parameters_min;
}

const RBParameters& RBParametrized::get_parameters_max() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameters_max" << std::endl;
    libmesh_error();
  }

  return parameters_max;
}

Real RBParametrized::get_parameter_min(const std::string& param_name) const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameter_min" << std::endl;
    libmesh_error();
  }

  return parameters_min.get_value(param_name);
}

Real RBParametrized::get_parameter_max(const std::string& param_name) const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::get_parameter_max" << std::endl;
    libmesh_error();
  }

  return parameters_max.get_value(param_name);
}

void RBParametrized::print_parameters() const
{
  if(!parameters_initialized)
  {
    libMesh::err << "Error: parameters not initialized in RBParametrized::print_current_parameters" << std::endl;
    libmesh_error();
  }

  get_parameters().print();
}

void RBParametrized::write_parameter_ranges_to_file(const std::string& file_name,
                                                    const bool write_binary_data)
{
  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = write_binary_data ? ENCODE : WRITE;

  // Write out the parameter ranges
  Xdr parameter_ranges_out(file_name, mode);
  unsigned int n_params = get_n_params();
  parameter_ranges_out << n_params;

  RBParameters::const_iterator it;
  RBParameters::const_iterator it_end;
  it = get_parameters_min().begin();
  it_end = get_parameters_min().end();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    Real param_value = it->second;

    parameter_ranges_out << param_name << param_value;
  }
  it     = get_parameters_max().begin();
  it_end = get_parameters_max().end();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    Real param_value = it->second;

    parameter_ranges_out << param_name << param_value;
  }
  parameter_ranges_out.close();
}

void RBParametrized::read_parameter_ranges_from_file(const std::string& file_name,
                                                     const bool read_binary_data)
{
  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // Read in the parameter ranges
  Xdr parameter_ranges_in(file_name, mode);
  unsigned int n_params;
  parameter_ranges_in >> n_params;
  RBParameters param_min;
  for(unsigned int i=0; i<n_params; i++)
  {
    std::string param_name;
    Real param_value;

    parameter_ranges_in >> param_name;
    parameter_ranges_in >> param_value;

    param_min.set_value(param_name, param_value);
  }
  RBParameters param_max;
  for(unsigned int i=0; i<n_params; i++)
  {
    std::string param_name;
    Real param_value;

    parameter_ranges_in >> param_name;
    parameter_ranges_in >> param_value;

    param_max.set_value(param_name, param_value);
  }
  parameter_ranges_in.close();

  initialize_parameters(param_min, param_max, param_min);
}

bool RBParametrized::valid_params(const RBParameters& params)
{
  if(params.n_parameters() != get_n_params())
  {
    libMesh::out << "Error: Number of parameters don't match" << std::endl;
    libmesh_error();
    return false;
  }
  else
  {
    bool valid = true;
    RBParameters::const_iterator it     = params.begin();
    RBParameters::const_iterator it_end = params.end();
    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      valid = valid && ( (get_parameter_min(param_name) <= params.get_value(param_name)) &&
                         (params.get_value(param_name) <= get_parameter_max(param_name)) );
    }

    if(!valid && verbose_mode)
    {
      libMesh::out << "Warning: parameter is outside parameter range" << std::endl;
    }

    return valid;
  }
}

} // namespace libMesh
