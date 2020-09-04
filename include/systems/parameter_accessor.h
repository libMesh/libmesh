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



#ifndef LIBMESH_PARAMETER_ACCESSOR_H
#define LIBMESH_PARAMETER_ACCESSOR_H


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/compare_types.h" // remove_const

// C++ includes
#include <memory>

namespace libMesh
{

// Forward declarations
template <typename T>
class ParameterProxy;

template <typename T>
class ConstParameterProxy;


/**
 * Accessor object allowing reading and modification of the
 * independent variables in a parameter sensitivity calculation.
 *
 * This is an abstract base class.  Derived objects may simply modify
 * the parameter value at some address in memory, or may call
 * arbitrary setter/getter functions.
 *
 * \author Roy Stogner
 * \date 2015
 * \brief Base class for reading/writing sensitivity parameters.
 */
template <typename T=Number>
class ParameterAccessor
{
public:
  /**
   * Virtual destructor - we'll be deleting subclasses from
   * pointers-to-ParameterAccessor
   */
  virtual ~ParameterAccessor() {}

  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const T & new_value) = 0;

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const = 0;

  /**
   * Proxy: for backward compatibility, we allow codes to treat a
   * ParameterAccessor as if it were a simple pointer-to-value.  We
   * can't safely allow "Number * n = parameter_vector[p]" to compile,
   * but we can allow "*parameter_vector[p] += deltap" to work.
   */
  ParameterProxy<T> operator* () { return ParameterProxy<T>(*this); }

  ConstParameterProxy<T> operator* () const { return ConstParameterProxy<T>(*this); }

  /**
   * \returns A new copy of the accessor.  The new copy should probably
   * be as shallow as possible, but should still access the same
   * parameter.
   */
  virtual std::unique_ptr<ParameterAccessor<T>> clone() const = 0;
};

template <typename T=Number>
class ParameterProxy
{
public:
  /**
   * Constructor: which parameter are we a proxy for?
   */
  ParameterProxy (ParameterAccessor<T> & accessor)
    : _accessor(accessor) {}

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator = (const T & new_value) { _accessor.set(new_value); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator = (const ParameterProxy<T> & new_value) { _accessor.set(new_value.get()); }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator = (const ConstParameterProxy<T> & new_value) { _accessor.set(new_value.get()); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator += (const T & value_increment) { _accessor.set(_accessor.get() + value_increment); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator -= (const T & value_decrement) { _accessor.set(_accessor.get() - value_decrement); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator *= (const T & value_multiplier) { _accessor.set(_accessor.get() * value_multiplier); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy & operator /= (const T & value_divisor) { _accessor.set(_accessor.get() / value_divisor); return *this; }

  /**
   * Getter: get the value of the parameter we access.
   */
  operator T () const { return _accessor.get(); }

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
  operator boost::multiprecision::backends::float128_backend () const { return _accessor.get().backend(); }
#endif

private:
  ParameterAccessor<T> & _accessor;
};


template <typename T=Number>
class ConstParameterProxy
{
public:
  /**
   * Constructor: which parameter are we a proxy for?
   */
  ConstParameterProxy (const ParameterAccessor<T> & accessor)
    : _accessor(accessor) {}

  /**
   * Getter: get the value of the parameter we access.
   */
  operator T () const { return _accessor.get(); }

  /**
   * Getter: get the value of the parameter we access.
   */
  T get() const { return _accessor.get(); }

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
  operator boost::multiprecision::backends::float128_backend () const { return _accessor.get().backend(); }
#endif

private:
  const ParameterAccessor<T> & _accessor;
};


} // namespace libMesh

#endif // LIBMESH_PARAMETER_ACCESSOR_H
