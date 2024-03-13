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

// libmesh includes
#include "libmesh/int_range.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/rb_parameters.h"
#include "libmesh/utility.h"

// C++ includes
#include <algorithm>
#include <sstream>

namespace libMesh
{

RBParameters::RBParameters() :
  _n_samples(1)
{
}

RBParameters::RBParameters(const std::map<std::string, Real> & parameter_map) :
  _n_samples(1)
{
  // Backwards compatible support for constructing an RBParameters
  // object from a map<string, Real>. We store a single entry in each
  // vector in the map.
  for (const auto & [key, val] : parameter_map)
    _parameters[key] = {{val}};
}

void RBParameters::clear()
{
  _n_samples = 1;
  _parameters.clear();
  _extra_parameters.clear();
}

bool RBParameters::has_value(const std::string & param_name) const
{
  return _parameters.count(param_name);
}

bool RBParameters::has_extra_value(const std::string & param_name) const
{
  return _extra_parameters.count(param_name);
}

Real RBParameters::get_value(const std::string & param_name) const
{
  // Simply return the [0]th entry of the vector if possible, otherwise error.
  libmesh_error_msg_if(this->n_samples() != 1,
    "Requesting value for parameter " << param_name << ", but parameter contains multiple samples.");
  return this->get_sample_value(param_name, /*sample_idx=*/0);
}

const RBParameter& RBParameters::get_vector_value(const std::string & param_name) const
{
  // Simply return the [0]th entry of the vector if possible, otherwise error.
  libmesh_error_msg_if(this->n_samples() != 1,
    "Requesting value for parameter " << param_name << ", but parameter contains multiple sample.");
  return this->get_sample_vector_value(param_name, /*sample_idx=*/0);
}

Real RBParameters::get_value(const std::string & param_name, const Real & default_val) const
{
  // Simply return the [0]th entry of the vector if possible, otherwise error.
  libmesh_error_msg_if(this->n_samples() != 1,
    "Requesting value for parameter " << param_name << ", but parameter contains multiple samples.");
  return this->get_sample_value(param_name, /*sample_idx=*/0, default_val);
}

const RBParameter & RBParameters::get_vector_value(const std::string & param_name, const RBParameter & default_val) const
{
  // Simply return the [0]th entry of the vector if possible, otherwise error.
  libmesh_error_msg_if(this->n_samples() != 1,
    "Requesting value for parameter " << param_name << ", but parameter contains multiple samples.");
  return this->get_sample_vector_value(param_name, /*sample_idx=*/0, default_val);
}

Real RBParameters::get_sample_value(const std::string & param_name, std::size_t sample_idx) const
{
  const auto & sample_vec = libmesh_map_find(_parameters, param_name);
  libmesh_error_msg_if(sample_idx >= sample_vec.size(), "Error getting value for parameter " << param_name);
  libmesh_error_msg_if(sample_vec[sample_idx].size() != 1,
    "Requesting Real value for parameter " << param_name << ", but parameter contains multiple values.");
  return sample_vec[sample_idx][0];
}

const RBParameter & RBParameters::get_sample_vector_value(const std::string & param_name, std::size_t sample_idx) const
{
  const auto & sample_vec = libmesh_map_find(_parameters, param_name);
  libmesh_error_msg_if(sample_idx >= sample_vec.size(), "Error getting value for parameter " << param_name);
  return sample_vec[sample_idx];
}

Real RBParameters::get_sample_value(const std::string & param_name, std::size_t sample_idx, const Real & default_val) const
{
  auto it = _parameters.find(param_name);
  if (it == _parameters.end() || sample_idx >= it->second.size())
    return default_val;
  libmesh_error_msg_if(it->second[sample_idx].size() != 1,
    "Requesting Real value for parameter " << param_name << ", but parameter contains multiple values.");
  return it->second[sample_idx][0];
}

const RBParameter & RBParameters::get_sample_vector_value(const std::string & param_name, std::size_t sample_idx, const RBParameter & default_val) const
{
  auto it = _parameters.find(param_name);
  if (it == _parameters.end() || sample_idx >= it->second.size())
    return default_val;
  return it->second[sample_idx];
}

void RBParameters::set_value(const std::string & param_name, Real value)
{
  // This version of set_value() does not take an index and is provided
  // for backwards compatibility. It creates a vector entry for the specified
  // param_name, overwriting any value(s) that were present.
  _parameters[param_name] = {{value}};
}

void RBParameters::set_value(const std::string & param_name, const RBParameter & value)
{
  // This version of set_value() does not take an index and is provided
  // for backwards compatibility. It creates a vector entry for the specified
  // param_name, overwriting any value(s) that were present.
  _parameters[param_name] = {value};
}

void
RBParameters::set_value_helper(std::map<std::string, std::vector<RBParameter>> & map,
                               const std::string & param_name,
                               const std::size_t index,
                               RBParameter value)
{
  // Get reference to vector of values for this parameter, creating it
  // if it does not already exist.
  auto & sample_vec = map[param_name];

  // If vec is already big enough, just set the value
  if (sample_vec.size() > index)
    sample_vec[index] = std::move(value);

  // Otherwise push_back() if the vec is just barely not big enough
  else if (sample_vec.size() == index)
    sample_vec.emplace_back(std::move(value));

  // Otherwise, allocate more space (padding with 0s) if vector is not
  // big enough to fit the user's requested index.
  else
    {
      RBParameter zero_parameter(value.size(), 0.0);
      sample_vec.resize(index+1, zero_parameter);
      sample_vec[index] = std::move(value);
    }
}

void RBParameters::set_value(const std::string & param_name, std::size_t index, Real value)
{
  this->set_value_helper(_parameters, param_name, index, {value});
}

void RBParameters::set_value(const std::string & param_name, std::size_t index, const RBParameter & value)
{
  this->set_value_helper(_parameters, param_name, index, value);
}

void RBParameters::set_extra_value(const std::string & param_name, std::size_t index, Real value)
{
  this->set_value_helper(_extra_parameters, param_name, index, {value});
}

void RBParameters::set_extra_value(const std::string & param_name, std::size_t index, const RBParameter & value)
{
  this->set_value_helper(_extra_parameters, param_name, index, value);
}

void RBParameters::push_back_value(const std::string & param_name, Real value)
{
  // Get reference to vector of values for this parameter, creating it
  // if it does not already exist, and push back the specified value.
  _parameters[param_name].push_back({value});
}

void RBParameters::push_back_value(const std::string & param_name, const RBParameter & value)
{
  // Get reference to vector of values for this parameter, creating it
  // if it does not already exist, and push back the specified value.
  _parameters[param_name].push_back(value);
}

void RBParameters::push_back_extra_value(const std::string & param_name, Real value)
{
  // Get reference to vector of values for this extra parameter, creating it
  // if it does not already exist, and push back the specified value.
  _extra_parameters[param_name].push_back({value});
}

void RBParameters::push_back_extra_value(const std::string & param_name, const RBParameter & value)
{
  // Get reference to vector of values for this extra parameter, creating it
  // if it does not already exist, and push back the specified value.
  _extra_parameters[param_name].push_back(value);
}

Real RBParameters::get_extra_value(const std::string & param_name) const
{
  // Same as get_value(param_name) but for the map of extra parameters
  const auto & sample_vec = libmesh_map_find(_extra_parameters, param_name);
  libmesh_error_msg_if(sample_vec.size() != 1,
    "Requesting value for extra parameter " << param_name << ", but parameter contains multiple samples.");
  libmesh_error_msg_if(sample_vec[0].size() != 1,
    "Requesting Real value for extra parameter " << param_name << ", but parameter contains multiple values.");
  return sample_vec[0][0];
}

const RBParameter & RBParameters::get_extra_vector_value(const std::string & param_name) const
{
  // Same as get_value(param_name) but for the map of extra parameters
  const auto & sample_vec = libmesh_map_find(_extra_parameters, param_name);
  libmesh_error_msg_if(sample_vec.size() == 0, "Error getting value for extra parameter " << param_name);
  return sample_vec[0];
}

Real RBParameters::get_extra_value(const std::string & param_name, const Real & default_val) const
{
  // same as get_value(param_name, default_val) but for the map of extra parameters
  auto it = _extra_parameters.find(param_name);
  if (it == _extra_parameters.end())
    return default_val;

  libmesh_error_msg_if(it->second.size() != 1,
                       "Requesting value for extra parameter " << param_name << ", but parameter contains multiple samples.");
  libmesh_error_msg_if(it->second[0].size() != 1,
    "Requesting Real value for extra parameter " << param_name << ", but parameter contains multiple values.");
  return it->second[0][0];
}

Real RBParameters::get_extra_sample_value(const std::string & param_name, std::size_t sample_idx) const
{
  const auto & sample_vec = libmesh_map_find(_extra_parameters, param_name);
  libmesh_error_msg_if(sample_idx >= sample_vec.size(), "Error getting value for extra parameter " << param_name);
  libmesh_error_msg_if(sample_vec[sample_idx].size() != 1,
    "Requesting Real value for extra parameter " << param_name << ", but parameter contains multiple values.");
  return sample_vec[sample_idx][0];
}

const RBParameter & RBParameters::get_extra_sample_vector_value(const std::string & param_name, std::size_t sample_idx) const
{
  const auto & sample_vec = libmesh_map_find(_extra_parameters, param_name);
  libmesh_error_msg_if(sample_idx >= sample_vec.size(), "Error getting value for extra parameter " << param_name);
  return sample_vec[sample_idx];
}

Real RBParameters::get_extra_sample_value(const std::string & param_name, std::size_t sample_idx, const Real & default_val) const
{
  // same as get_sample_value(param_name, index, default_val) but for the map of extra parameters
  auto it = _extra_parameters.find(param_name);
  if (it==_extra_parameters.end() || sample_idx >= it->second.size())
    return default_val;
  libmesh_error_msg_if(it->second[sample_idx].size() != 1,
    "Requesting Real value for extra parameter " << param_name << ", but parameter contains multiple values.");
  return it->second[sample_idx][0];
}

const RBParameter & RBParameters::get_extra_sample_vector_value(
    const std::string &param_name,
    std::size_t sample_idx,
    const RBParameter &default_val) const
{
  auto it = _extra_parameters.find(param_name);
  if (it == _extra_parameters.end() || sample_idx >= it->second.size())
    return default_val;
  return it->second[sample_idx];
}

void RBParameters::set_extra_value(const std::string & param_name, Real value)
{
  // This version of set_extra_value() does not take an index and is provided
  // for backwards compatibility. It creates a vector entry for the specified
  // param_name, overwriting any value(s) that were present.
  _extra_parameters[param_name] = {{value}};
}

void RBParameters::set_extra_value(const std::string & param_name, const RBParameter & value)
{
  // This version of set_extra_value() does not take an index and is provided
  // for backwards compatibility. It creates a vector entry for the specified
  // param_name, overwriting any value(s) that were present.
  _extra_parameters[param_name] = {value};
}

unsigned int RBParameters::n_parameters() const
{
  return cast_int<unsigned int>(_parameters.size());
}

void RBParameters::set_n_samples(unsigned int n_samples)
{
  _n_samples = n_samples;
}

unsigned int RBParameters::n_samples() const
{
  // Quick return if there are no parameters
  if (_parameters.empty())
    return _n_samples;

  // If _parameters is not empty, we can check the number of samples in the first param
  auto size_first = _parameters.begin()->second.size();

#ifdef DEBUG
  // In debug mode, verify that all parameters have the same number of samples
  for (const auto & pr : _parameters)
    libmesh_assert_msg(pr.second.size() == size_first, "All parameters must have the same number of samples.");
#endif

  // If we made it here in DEBUG mode, then all parameters were
  // verified to have the same number of samples.
  return size_first;
}

std::set<std::string> RBParameters::get_parameter_names() const
{
  std::set<std::string> param_names;
  for (const auto & pr : _parameters)
    param_names.insert(pr.first);
  return param_names;
}

std::set<std::string> RBParameters::get_extra_parameter_names() const
{
  std::set<std::string> param_names;
  for (const auto & pr : _extra_parameters)
    param_names.insert(pr.first);
  return param_names;
}

void RBParameters::erase_parameter(const std::string & param_name)
{
  _parameters.erase(param_name);
}

void RBParameters::erase_extra_parameter(const std::string & param_name)
{
  _extra_parameters.erase(param_name);
}

std::map<std::string,std::vector<RBParameter>>::const_iterator RBParameters::begin() const
{
  return _parameters.cbegin();
}

std::map<std::string,std::vector<RBParameter>>::const_iterator RBParameters::end() const
{
  return _parameters.cend();
}

std::map<std::string,std::vector<RBParameter>>::const_iterator RBParameters::extra_begin() const
{
  return _extra_parameters.cbegin();
}

std::map<std::string,std::vector<RBParameter>>::const_iterator RBParameters::extra_end() const
{
  return _extra_parameters.cend();
}

RBParameters::const_iterator RBParameters::begin_serialized() const
{
  return {_parameters.cbegin(), 0, 0};
}

RBParameters::const_iterator RBParameters::end_serialized() const
{
  // Note: the index 0 is irrelevant here since _parameters.end() does
  // not refer to a valid vector entry in the map.
  return {_parameters.cend(), 0, 0};
}

RBParameters::const_iterator RBParameters::begin_serialized_extra() const
{
  return {_extra_parameters.cbegin(), 0, 0};
}

RBParameters::const_iterator RBParameters::end_serialized_extra() const
{
  // Note: the index 0 is irrelevant here since _parameters.end() does
  // not refer to a valid vector entry in the map.
  return {_extra_parameters.cend(), 0, 0};
}

bool RBParameters::operator==(const RBParameters & rhs) const
{
  return (this->_parameters == rhs._parameters &&
          this->_extra_parameters == rhs._extra_parameters);
}

bool RBParameters::operator!=(const RBParameters & rhs) const
{
  return !(*this == rhs);
}

RBParameters & RBParameters::operator+= (const RBParameters & rhs)
{
  libmesh_error_msg_if(this->n_samples() != rhs.n_samples(),
                       "Can only append RBParameters objects with matching numbers of samples.");

  // Overwrite or add each (key, vec) pair in rhs to *this.
  for (const auto & [key, vec] : rhs._parameters)
    _parameters[key] = vec;
  for (const auto & [key, vec] : rhs._extra_parameters)
    _extra_parameters[key] = vec;

  return *this;
}

std::string RBParameters::get_string(unsigned precision, int max_values) const
{
  std::stringstream param_stringstream;
  param_stringstream << std::setprecision(static_cast<int>(precision)) << std::scientific;

  for (const auto & [param_name, sample_vec] : _parameters)
    {
      // Write the param name, followed by a comma-separated list of the sample/vector values.
      param_stringstream << param_name << ": ";
      std::string separator = "";
      for (const auto & value_vec : sample_vec)
        {
          param_stringstream << separator;
          if (value_vec.size() == 1)
            param_stringstream << value_vec[0];
          else
            {
              param_stringstream << "[ ";
              for (const auto val_idx : index_range(value_vec))
                {
                  if (max_values < 0 || val_idx < static_cast<unsigned>(max_values))
                    param_stringstream << value_vec[val_idx] << " ";
                  else
                    {
                      param_stringstream << "... ";
                      break;
                    }
                }
              param_stringstream << "]";
            }
          separator = ", ";
        }
      param_stringstream << std::endl;
    }

  return param_stringstream.str();
}

void RBParameters::print(unsigned precision, int max_values) const
{
  libMesh::out << get_string(precision, max_values);
}

}
