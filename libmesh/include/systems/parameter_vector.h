

// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef __parameter_vector_h__
#define __parameter_vector_h__


// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
#include "libmesh_common.h"


  /**
   * Data structure for specifying which Parameters should be
   * independent variables in a parameter sensitivity calculation.
   */
class ParameterVector
{
public: 
  /**
   * Default constructor: "no parameters"
   */
  ParameterVector() {}

  /**
   * Constructor-from-vector-of-Number*: each points to a parameter
   */
  ParameterVector(const std::vector<Number *> &params) : _params(params) {}

  /**
   * Resets to "no parameters"
   */
  void clear() { _params.clear(); }

  /**
   * Returns the number of parameters to be used
   */
  unsigned int size() const { return _params.size(); }

  /**
   * Sets the number of parameters to be used
   */
  void resize(unsigned int s) { _params.resize(s); }

  /**
   * Returns a pointer to a parameter value
   */
  Number * operator[](unsigned int i) const;

  /**
   * Returns a reference to a pointer to a parameter value,
   * suitable for repointing it to a different address.
   */
  Number *& operator[](unsigned int i);

private: 
  std::vector<Number *> _params;
};



// ------------------------------------------------------------
// ParameterVector inline methods



inline
Number* ParameterVector::operator[] (unsigned int i) const
{
  libmesh_assert(_params.size() > i);

  return _params[i];
}



inline
Number*& ParameterVector::operator[] (unsigned int i)
{
  libmesh_assert(_params.size() > i);

  return _params[i];
}

#endif // #define __parameter_vector_h__
