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
#include <fstream>
#include <algorithm> // std::min_element

namespace libMesh
{

// ------------------------------------------------------------
// RBParameters implementation

RBParametrized::RBParametrized()
  :
  verbose_mode(false),
  parameters_initialized(false)
{
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

void RBParametrized::initialize_parameters(const RBParameters & mu_min_in,
                                           const RBParameters & mu_max_in,
                                           const std::map< std::string, std::vector<Real> > & discrete_parameter_values)
{
  // Check that the min/max vectors are valid
  {
    const std::string err_string = "Error: Invalid mu_min/mu_max in RBParameters constructor.";
    bool valid_min_max = (mu_min_in.n_parameters() == mu_max_in.n_parameters());
    if(!valid_min_max)
      libmesh_error_msg(err_string);

    else
      {
        RBParameters::const_iterator it     = mu_min_in.begin();
        RBParameters::const_iterator it_end = mu_min_in.end();
        for( ; it != it_end; ++it)
          {
            std::string param_name = it->first;
            if(mu_min_in.get_value(param_name) > mu_max_in.get_value(param_name))
              libmesh_error_msg(err_string);
          }
      }
  }

  parameters_min = mu_min_in;
  parameters_max = mu_max_in;

  // Add in min/max values due to the discrete parameters
  {
    std::map< std::string, std::vector<Real> >::const_iterator it =
      discrete_parameter_values.begin();
    const std::map< std::string, std::vector<Real> >::const_iterator it_end =
      discrete_parameter_values.end();
    for( ; it != it_end; ++it)
      {
        if(it->second.empty())
          libmesh_error_msg("Error: List of discrete parameters for " << it->first << " is empty.");

        Real min_val = *std::min_element(it->second.begin(), it->second.end());
        Real max_val = *std::max_element(it->second.begin(), it->second.end());

        libmesh_assert_less_equal(min_val,max_val);

        parameters_min.set_value(it->first, min_val);
        parameters_max.set_value(it->first, max_val);
      }

    _discrete_parameter_values = discrete_parameter_values;
  }

  parameters_initialized = true;

  // Initialize the current parameters to parameters_min
  set_parameters(parameters_min);
}

void RBParametrized::initialize_parameters(const RBParametrized & rb_parametrized)
{
  initialize_parameters(rb_parametrized.get_parameters_min(),
                        rb_parametrized.get_parameters_max(),
                        rb_parametrized.get_discrete_parameter_values());
}

unsigned int RBParametrized::get_n_params() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_n_params");

  libmesh_assert_equal_to ( parameters_min.n_parameters(), parameters_max.n_parameters() );

  return parameters_min.n_parameters();
}

unsigned int RBParametrized::get_n_continuous_params() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_n_continuous_params");

  libmesh_assert(get_n_params() >= get_n_discrete_params());

  return static_cast<unsigned int>(get_n_params() - get_n_discrete_params());
}

unsigned int RBParametrized::get_n_discrete_params() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_n_discrete_params");

  return cast_int<unsigned int>
    (get_discrete_parameter_values().size());
}

std::set<std::string> RBParametrized::get_parameter_names() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_parameter_names");

  std::set<std::string> parameter_names;
  parameters_min.get_parameter_names(parameter_names);

  return parameter_names;
}

void RBParametrized::set_parameters(const RBParameters & params)
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::set_current_parameters");

  valid_params(params); // Terminates if params has the wrong number of parameters

  // Make a copy of params (default assignment operator just does memberwise copy, which is sufficient here)
  this->parameters = params;
}

const RBParameters & RBParametrized::get_parameters() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_current_parameters");

  return parameters;
}

const RBParameters & RBParametrized::get_parameters_min() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_parameters_min");

  return parameters_min;
}

const RBParameters & RBParametrized::get_parameters_max() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_parameters_max");

  return parameters_max;
}

Real RBParametrized::get_parameter_min(const std::string & param_name) const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_parameter_min");

  return parameters_min.get_value(param_name);
}

Real RBParametrized::get_parameter_max(const std::string & param_name) const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_parameter_max");

  return parameters_max.get_value(param_name);
}

void RBParametrized::print_parameters() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::print_current_parameters");

  get_parameters().print();
}

void RBParametrized::write_parameter_data_to_files(const std::string & continuous_param_file_name,
                                                   const std::string & discrete_param_file_name,
                                                   const bool write_binary_data)
{
  write_parameter_ranges_to_file(continuous_param_file_name, write_binary_data);
  write_discrete_parameter_values_to_file(discrete_param_file_name, write_binary_data);
}

void RBParametrized::write_parameter_ranges_to_file(const std::string & file_name,
                                                    const bool write_binary_data)
{
  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = write_binary_data ? ENCODE : WRITE;

  // Write out the parameter ranges
  Xdr parameter_ranges_out(file_name, mode);
  unsigned int n_continuous_params = get_n_continuous_params();
  parameter_ranges_out << n_continuous_params;

  RBParameters::const_iterator it;
  RBParameters::const_iterator it_end;
  it = get_parameters_min().begin();
  it_end = get_parameters_min().end();
  for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      if(!is_discrete_parameter(param_name))
        {
          Real param_value = it->second;
          parameter_ranges_out << param_name << param_value;
        }
    }
  it     = get_parameters_max().begin();
  it_end = get_parameters_max().end();
  for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      if(!is_discrete_parameter(param_name))
        {
          Real param_value = it->second;
          parameter_ranges_out << param_name << param_value;
        }
    }
  parameter_ranges_out.close();
}

void RBParametrized::write_discrete_parameter_values_to_file(const std::string & file_name,
                                                             const bool write_binary_data)
{
  // write out the discrete parameters, if we have any
  if(get_n_discrete_params() > 0)
    {
      // The writing mode: ENCODE for binary, WRITE for ASCII
      XdrMODE mode = write_binary_data ? ENCODE : WRITE;

      Xdr discrete_parameters_out(file_name, mode);
      unsigned int n_discrete_params = get_n_discrete_params();
      discrete_parameters_out << n_discrete_params;

      std::map< std::string, std::vector<Real> >::const_iterator discrete_it =
        get_discrete_parameter_values().begin();
      const std::map< std::string, std::vector<Real> >::const_iterator discrete_it_end =
        get_discrete_parameter_values().end();
      for( ; discrete_it != discrete_it_end; ++discrete_it)
        {
          std::string param_name = discrete_it->first;
          unsigned int n_discrete_values = cast_int<unsigned int>
            (discrete_it->second.size());
          discrete_parameters_out << param_name << n_discrete_values;

          for(unsigned int i=0; i<n_discrete_values; i++)
            {
              Real discrete_value = discrete_it->second[i];
              discrete_parameters_out << discrete_value;
            }
        }
    }
}

void RBParametrized::read_parameter_data_from_files(const std::string & continuous_param_file_name,
                                                    const std::string & discrete_param_file_name,
                                                    const bool read_binary_data)
{
  RBParameters param_min;
  RBParameters param_max;
  read_parameter_ranges_from_file(continuous_param_file_name,
                                  read_binary_data,
                                  param_min,
                                  param_max);

  std::map< std::string, std::vector<Real> > discrete_parameter_values_in;
  read_discrete_parameter_values_from_file(discrete_param_file_name,
                                           read_binary_data,
                                           discrete_parameter_values_in);

  initialize_parameters(param_min, param_max, discrete_parameter_values_in);
}

void RBParametrized::read_parameter_ranges_from_file(const std::string & file_name,
                                                     const bool read_binary_data,
                                                     RBParameters & param_min,
                                                     RBParameters & param_max)
{
  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // Read in the parameter ranges
  Xdr parameter_ranges_in(file_name, mode);
  unsigned int n_continuous_params;
  parameter_ranges_in >> n_continuous_params;

  for(unsigned int i=0; i<n_continuous_params; i++)
    {
      std::string param_name;
      Real param_value;

      parameter_ranges_in >> param_name;
      parameter_ranges_in >> param_value;

      param_min.set_value(param_name, param_value);
    }
  for(unsigned int i=0; i<n_continuous_params; i++)
    {
      std::string param_name;
      Real param_value;

      parameter_ranges_in >> param_name;
      parameter_ranges_in >> param_value;

      param_max.set_value(param_name, param_value);
    }

  parameter_ranges_in.close();
}

void RBParametrized::read_discrete_parameter_values_from_file(const std::string & file_name,
                                                              const bool read_binary_data,
                                                              std::map< std::string, std::vector<Real> > & discrete_parameter_values)
{
  // read in the discrete parameters, if we have any
  std::ifstream check_if_file_exists(file_name.c_str());
  if(check_if_file_exists.good())
    {
      // The reading mode: DECODE for binary, READ for ASCII
      XdrMODE mode = read_binary_data ? DECODE : READ;

      // Read in the parameter ranges
      Xdr discrete_parameter_values_in(file_name, mode);
      unsigned int n_discrete_params;
      discrete_parameter_values_in >> n_discrete_params;

      for(unsigned int i=0; i<n_discrete_params; i++)
        {
          std::string param_name;
          discrete_parameter_values_in >> param_name;

          unsigned int n_discrete_values;
          discrete_parameter_values_in >> n_discrete_values;

          std::vector<Real> discrete_values(n_discrete_values);
          for (std::size_t j=0; j<discrete_values.size(); j++)
            {
              Real param_value;
              discrete_parameter_values_in >> param_value;

              discrete_values[j] = param_value;
            }

          discrete_parameter_values[param_name] = discrete_values;
        }
    }
}

bool RBParametrized::is_discrete_parameter(const std::string & mu_name) const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::is_discrete_parameter");

  return (_discrete_parameter_values.find(mu_name) != _discrete_parameter_values.end());
}

const std::map< std::string, std::vector<Real> > & RBParametrized::get_discrete_parameter_values() const
{
  if(!parameters_initialized)
    libmesh_error_msg("Error: parameters not initialized in RBParametrized::get_discrete_parameter_values");

  return _discrete_parameter_values;
}

void RBParametrized::print_discrete_parameter_values() const
{
  std::map< std::string, std::vector<Real> >::const_iterator it =
    get_discrete_parameter_values().begin();
  const std::map< std::string, std::vector<Real> >::const_iterator it_end =
    get_discrete_parameter_values().end();

  for( ; it != it_end; ++it)
    {
      libMesh::out << "Discrete parameter " << it->first << ", values: ";

      std::vector<Real> values = it->second;
      for (std::size_t i=0; i<values.size(); i++)
        libMesh::out << values[i] << " ";
      libMesh::out << std::endl;
    }
}

bool RBParametrized::valid_params(const RBParameters & params)
{
  if(params.n_parameters() != get_n_params())
    libmesh_error_msg("Error: Number of parameters don't match");

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

          if(is_discrete_parameter(param_name))
            {
              // make sure params.get_value(param_name) is sufficiently close
              // to one of the discrete parameter values
              valid = valid && is_value_in_list(params.get_value(param_name),
                                                get_discrete_parameter_values().find(param_name)->second,
                                                TOLERANCE);
            }
        }

      if(!valid && verbose_mode)
        {
          libMesh::out << "Warning: parameter is outside parameter range" << std::endl;
        }

      return valid;
    }
}

Real RBParametrized::get_closest_value(Real value, const std::vector<Real> & list_of_values)
{
  if(list_of_values.empty())
    libmesh_error_msg("Error: list_of_values is empty.");

  std::vector<Real>::const_iterator it = list_of_values.begin();
  std::vector<Real>::const_iterator it_end = list_of_values.end();

  Real min_distance = std::numeric_limits<Real>::max();
  Real closest_val = 0.;
  for( ; it != it_end; ++it)
    {
      Real distance = std::abs(value - *it);
      if( distance < min_distance )
        {
          min_distance = distance;
          closest_val = *it;
        }
    }

  return closest_val;
}

bool RBParametrized::is_value_in_list(Real value, const std::vector<Real> & list_of_values, Real tol)
{
  Real closest_value = get_closest_value(value, list_of_values);

  // Check if relative tolerance is satisfied
  Real rel_error = std::abs(value - closest_value) / std::abs(value);
  if( rel_error <= tol )
    {
      return true;
    }

  // If relative tolerance isn't satisfied, we should still check an absolute
  // error, since relative tolerance can be misleading if value is close to zero
  Real abs_error = std::abs(value - closest_value);
  return (abs_error <= tol);
}

} // namespace libMesh
