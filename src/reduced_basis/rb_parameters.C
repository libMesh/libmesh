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

// C++ includes
#include <sstream>

// libmesh includes
#include "libmesh/rb_parameters.h"
#include "libmesh/utility.h"

namespace libMesh
{

RBParameters::RBParameters(const std::map<std::string, Real> & parameter_map)
{
  _parameters = parameter_map;
}

void RBParameters::clear()
{
  _parameters.clear();
  _extra_parameters.clear();
}

const std::map<std::string, Real> & RBParameters::get_parameters_map() const
{
  return _parameters;
}

const std::map<std::string, Real> & RBParameters::get_extra_parameters_map() const
{
  return _extra_parameters;
}

Real RBParameters::get_value(const std::string & param_name) const
{
  // find the parameter value, throwing an error if it doesn't exist.
  return libmesh_map_find(_parameters, param_name);
}

void RBParameters::set_value(const std::string & param_name, Real value)
{
  _parameters[param_name] = value;
}

Real RBParameters::get_extra_value(const std::string & param_name) const
{
  // find the parameter value, throwing an error if it doesn't exist.
  return libmesh_map_find(_extra_parameters, param_name);
}

void RBParameters::set_extra_value(const std::string & param_name, Real value)
{
  _extra_parameters[param_name] = value;
}

unsigned int RBParameters::n_parameters() const
{
  return cast_int<unsigned int>
    (_parameters.size());
}

void RBParameters::get_parameter_names(std::set<std::string> & param_names) const
{
  libmesh_deprecated();

  param_names.clear();
  for (const auto & pr : _parameters)
    param_names.insert(pr.first);
}

void RBParameters::get_extra_parameter_names(std::set<std::string> & param_names) const
{
  libmesh_deprecated();

  param_names.clear();
  for (const auto & pr : _extra_parameters)
    param_names.insert(pr.first);
}

void RBParameters::erase_parameter(const std::string & param_name)
{
  _parameters.erase(param_name);
}

void RBParameters::erase_extra_parameter(const std::string & param_name)
{
  _extra_parameters.erase(param_name);
}

RBParameters::const_iterator RBParameters::begin() const
{
  return _parameters.begin();
}

RBParameters::const_iterator RBParameters::end() const
{
  return _parameters.end();
}

RBParameters::const_iterator RBParameters::extra_begin() const
{
  return _extra_parameters.begin();
}

RBParameters::const_iterator RBParameters::extra_end() const
{
  return _extra_parameters.end();
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

std::string RBParameters::get_string(unsigned int precision) const
{
  std::stringstream param_stringstream;
  param_stringstream.precision(precision);

  const_iterator it     = _parameters.begin();
  const_iterator it_end = _parameters.end();
  for ( ; it != it_end; ++it)
    {
      param_stringstream << it->first << ": " << std::scientific <<  it->second << std::endl;
    }
  return param_stringstream.str();
}

void RBParameters::print() const
{
  libMesh::out << get_string() << std::endl;
}

}
