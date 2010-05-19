

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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "parameter_vector.h"

// ------------------------------------------------------------
// ParameterVector implementation



void ParameterVector::deep_copy(ParameterVector &target) const
{
  const unsigned int Np = this->_params.size();
  target._params.resize(Np);
  target._my_data.resize(Np);
  for (unsigned int i=0; i != Np; ++i)
    {
      target._params[i] = &target._my_data[i];
      target._my_data[i] = *(this->_params[i]);
    }
}



void ParameterVector::shallow_copy(ParameterVector &target) const
{
  target._my_data.clear();
  target._params = this->_params;
}



void ParameterVector::value_copy(const ParameterVector &target) const
{
  const unsigned int Np = this->_params.size();
  libmesh_assert(target._params.size() == Np);

  for (unsigned int i=0; i != Np; ++i)
    *(this->_params[i]) = *(target._params[i]);
}



void ParameterVector::deep_resize(unsigned int s)
{
  this->_params.resize(s);
  this->_my_data.resize(s);
  for (unsigned int i=0; i != s; ++i)
    this->_params[i] = &this->_my_data[i];
}



ParameterVector& ParameterVector::operator *= (const Number a)
{
  const unsigned int Np = this->_params.size();
  for (unsigned int i=0; i != Np; ++i)
    *(this->_params[i]) *= a;
  return *this;
}



const ParameterVector& ParameterVector::operator += (const ParameterVector& a) const
{
  const unsigned int Np = this->_params.size();
  libmesh_assert(a._params.size() == Np);
  for (unsigned int i=0; i != Np; ++i)
    *(this->_params[i]) += *(a._params[i]);
  return *this;
}


ParameterVector& ParameterVector::operator += (const ParameterVector& a)
{
  (*this) += a;
  return *this;
}
