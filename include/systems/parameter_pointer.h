// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_accessor.h"

namespace libMesh
{

/**
 * Accessor object allowing reading and modification of the
 * independent variables in a parameter sensitivity calculation.
 *
 * This is the "default" ParameterAccessor subclass: it simply stores
 * a user-provided pointer to the parameter, and modifies the value at
 * that location in memory.
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
  virtual void set (const T & new_value) libmesh_override
  { libmesh_assert(_ptr); *_ptr = new_value; }

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const libmesh_override
  { libmesh_assert(_ptr); return *_ptr; }

  /**
   * Reseater: change the location of the parameter we access.
   * This is included for backward compatibility, but is deprecated.
   */
  virtual ParameterAccessor<T> &
  operator= (T * new_ptr) libmesh_override
  {
    libmesh_deprecated();
    _ptr = new_ptr;
    return *this;
  }

  /**
   * Returns a new copy of the accessor.
   */
  virtual UniquePtr<ParameterAccessor<T> > clone() const libmesh_override
  {
    return UniquePtr<ParameterAccessor<T> >(new ParameterPointer<T>(_ptr));
  }

private:
  T * _ptr;
};

} // namespace libMesh

#endif // LIBMESH_PARAMETER_POINTER_H
