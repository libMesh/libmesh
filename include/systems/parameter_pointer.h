// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARAMETER_POINTER_H
#define LIBMESH_PARAMETER_POINTER_H


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_accessor.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

/**
 * Accessor object allowing reading and modification of the
 * independent variables in a parameter sensitivity calculation.
 *
 * This is the "default" ParameterAccessor subclass: it simply stores
 * a user-provided pointer to the parameter, and modifies the value at
 * that location in memory.
 *
 * \author Roy Stogner
 * \date 2015
 * \brief Stores/modifies a user-provided pointer to a parameter.
 */
template <typename T=Number>
class ParameterPointer : public ParameterAccessor<T>
{
public:
  /**
   * Constructor: take the raw pointer to the parameter
   */
  ParameterPointer(T * param_ptr) : _ptr(param_ptr) {}
  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const T & new_value) override
  { libmesh_assert(_ptr); *_ptr = new_value; }

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const override
  { libmesh_assert(_ptr); return *_ptr; }

  /**
   * \returns A new copy of the accessor.
   */
  virtual std::unique_ptr<ParameterAccessor<T>> clone() const override
  {
    return libmesh_make_unique<ParameterPointer<T>>(_ptr);
  }

private:
  T * _ptr;
};

} // namespace libMesh

#endif // LIBMESH_PARAMETER_POINTER_H
