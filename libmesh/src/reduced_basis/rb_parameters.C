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
#include "rb_parameters.h"

namespace libMesh
{

RBParameters::RBParameters()
{
}

void RBParameters::clear()
{
  _parameters.clear();
}

Real RBParameters::get_value(const std::string& param_name) const
{
  // find the parameter value
  const_iterator it = _parameters.find(param_name);
  
  // throw and error if the parameter doesn't exist
  if( it == _parameters.end() )
  {
    libMesh::err << "Error: parameter " << param_name << " does not exist in RBParameters object."
                 << std::endl;
    libmesh_error();
  }
  
  return it->second;
}

void RBParameters::set_value(const std::string& param_name, Real value)
{
  _parameters[param_name] = value;
}

unsigned int RBParameters::n_parameters() const
{
  return _parameters.size();
}

void RBParameters::get_parameter_names(std::set<std::string>& param_names) const
{
  param_names.clear();
  
  const_iterator it     = _parameters.begin();
  const_iterator it_end = _parameters.end();
  for( ; it != it_end; ++it)
  {
    param_names.insert( it->first );
  }
}

RBParameters::const_iterator RBParameters::begin() const
{
  return _parameters.begin();
}

RBParameters::const_iterator RBParameters::end() const
{
  return _parameters.end();
}

bool RBParameters::operator==(const RBParameters& rhs) const
{
  return this->_parameters == rhs._parameters;
}

bool RBParameters::operator!=(const RBParameters& rhs) const
{
  return !(*this == rhs);
}

std::string RBParameters::get_string() const
{
  std::stringstream param_stringstream;
  
  const_iterator it     = _parameters.begin();
  const_iterator it_end = _parameters.end();
  for( ; it != it_end; ++it)
  {
    param_stringstream << it->first << ": " << it->second << std::endl;
  }
  return param_stringstream.str();
}

void RBParameters::print() const
{
  libMesh::out << get_string() << std::endl;
}

}
