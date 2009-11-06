

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
   * Deep copy constructor: the \p target will now own new copies of
   * all the parameter values I'm pointing to
   */
  void deep_copy(ParameterVector &target) const;

  /**
   * Shallow copy constructor: the \p target will now point to all the
   * parameter values I'm pointing to
   */
  void shallow_copy(ParameterVector &target) const;

  /**
   * Value copy method: the \p target, which should already have as
   * many parameters as I do, will now have those parameters set to my
   * values.
   */
  void value_copy(const ParameterVector &target) const;

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

  /**
   * Multiplication operator; acts individually on each parameter.
   */
  ParameterVector& operator *= (const Number a);

  /**
   * Addition operator.  The parameter vector to be added in must
   * have the same number of values.
   */
  ParameterVector& operator += (const ParameterVector& a);

  /**
   * Addition operator.  The parameter vector to be added in must
   * have the same number of values.
   */
  const ParameterVector& operator += (const ParameterVector& a) const;

private: 
  /**
   * Pointers to parameters which may exist elsewhere
   */
  std::vector<Number *> _params;

  /**
   * Parameters which I own; e.g. as the result of a deep copy
   */
  std::vector<Number> _my_data;
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
