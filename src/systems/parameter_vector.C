// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes
#include "libmesh/parameter_vector.h"

#include "libmesh/int_range.h"
#include "libmesh/parameter_pointer.h"

namespace libMesh
{

ParameterVector::ParameterVector(const std::vector<Number *> &params)
#ifndef NDEBUG
  : _is_shallow_copy(false)
#endif
{
  _params.reserve(params.size());

  for (auto p : params)
    _params.push_back(std::make_unique<ParameterPointer<Number>>(p));
}



void ParameterVector::deep_copy(ParameterVector & target) const
{
  const std::size_t Np = this->_params.size();
  target.clear();
  target._params.resize(Np);
  target._my_data.resize(Np);
#ifndef NDEBUG
  target._is_shallow_copy = false;
#endif
  for (std::size_t i=0; i != Np; ++i)
    {
      target._params[i] =
        std::make_unique<ParameterPointer<Number>>
          (&target._my_data[i]);
      target._my_data[i] = *(*this)[i];
    }
}



void ParameterVector::shallow_copy(ParameterVector & target) const
{
  target._my_data.clear();
  target._params.resize(this->_params.size());
  for (auto i : index_range(this->_params))
    target._params[i] = this->_params[i]->clone();
#ifndef NDEBUG
  target._is_shallow_copy = true;
#endif
}



void ParameterVector::value_copy(ParameterVector & target) const
{
  const std::size_t Np = this->_params.size();
  libmesh_assert_equal_to (target._params.size(), Np);

  for (std::size_t i=0; i != Np; ++i)
    *target[i] = *(*this)[i];
}



void ParameterVector::resize(std::size_t s)
{
  libmesh_assert(!_is_shallow_copy);

  this->_params.resize(s);

  // We used to make nullptr ParameterPointers here, but that was just
  // to keep our destructor happy; users couldn't reseat them so
  // shouldn't have code depending on them
}



void ParameterVector::deep_resize(std::size_t s)
{
  libmesh_assert(!_is_shallow_copy);

  this->_params.resize(s);
  this->_my_data.resize(s);
  for (std::size_t i=0; i != s; ++i)
    this->_params[i] =
      std::make_unique<ParameterPointer<Number>>(&this->_my_data[i]);
}



ParameterVector & ParameterVector::operator *= (const Number a)
{
  const std::size_t Np = this->_params.size();
  for (std::size_t i=0; i != Np; ++i)
    *(*this)[i] *= a;
  return *this;
}



ParameterVector & ParameterVector::operator += (const ParameterVector & a)
{
  const std::size_t Np = this->_params.size();
  libmesh_assert_equal_to (a._params.size(), Np);
  for (std::size_t i=0; i != Np; ++i)
    *(*this)[i] += *a[i];
  return *this;
}


} // namespace libMesh
