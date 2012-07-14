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

#ifndef __rb_parameters_h__
#define __rb_parameters_h__

// libMesh includes
#include "libmesh_common.h"

// C++ includes
#include <string>
#include <map>
#include <set>

namespace libMesh
{
    
/**
 * This class is part of the rbOOmit framework.
 *
 * This class defines a set of parameters index by strings.
 *
 * @author David J. Knezevic, 2012
 */
// ------------------------------------------------------------
// RBParameters class definition
class RBParameters
{
public:

  /**
   * Constructor.
   */
  RBParameters();

  // Define a constant iterator for this class
  typedef std::map<std::string, Real>::const_iterator const_iterator;

  /**
   * Clear this object.
   */
  void clear();
  
  /**
   * Get the value of the specific parameter.
   */
  Real get_value(const std::string& param_name) const;
  
  /**
   * Set the value of the specified parameter. If param_name
   * doesn't already exist, it is added to the RBParameters object.
   */
  void set_value(const std::string& param_name, Real value);
  
  /**
   * Get the number of parameters that have been added.
   */
  unsigned int n_parameters() const;
  
  /**
   * Fill \p param_names with the names of the parameters.
   */
  void get_parameter_names(std::set<std::string>& param_names) const;
  
  /**
   * Get a constant iterator to beginning of this RBParameters object.
   */
  const_iterator begin() const;
  
  /**
   * Get a constant iterator to the end of this RBParameters object.
   */
  const_iterator end() const;
  
  /**
   * Two RBParameters are equal if they have the same _parameters map.
   */
  bool operator== (const RBParameters& rhs) const;

  /**
   * @return !(*this == rhs).
   */
  bool operator!= (const RBParameters& node) const;
  
  /**
   * Print the parameters.
   */
  void print() const;

private:

  /**
   * The map that stores the actual parameters, indexed by names.
   */
  std::map<std::string, Real> _parameters;

};

} // namespace libMesh


#endif
