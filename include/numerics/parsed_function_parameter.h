// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARSED_FUNCTION_PARAMETER_H
#define LIBMESH_PARSED_FUNCTION_PARAMETER_H


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_accessor.h"
#include "libmesh/parsed_function.h"

namespace libMesh
{

/**
 * Accessor object allowing reading and modification of the
 * independent variables in a parameter sensitivity calculation.
 *
 * This ParameterAccessor subclass is specific to ParsedFunction
 * objects: it stores a pointer to the ParsedFunction and a string
 * describing the parameter (an inline variable) name to be accessed.
 */
template <typename T=Number>
class ParsedFunctionParameter : public ParameterAccessor<T>
{
public:
  /**
   * Constructor: take the function to be modified and the name of the
   * inline variable within it which represents our parameter.
   *
   * The restrictions of get_inline_value() and set_inline_value()
   * in ParsedFunction apply to this interface as well.
   *
   * Note that *only* the function referred to here is changed by
   * set() - any clones of the function which precede the set()
   * remain at their previous values.
   */
  ParsedFunctionParameter(ParsedFunction<T> & func_ref,
                          const std::string & param_name) :
    _func(func_ref), _name(param_name) {}

  /**
   * A simple reseater won't work with a parsed function
   */
  virtual ParameterAccessor<T> &
  operator= (T * /* new_ptr */) { libmesh_error(); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const T & new_value) {
    _func.set_inline_value(_name, new_value);
  }

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const {
    _current_val = _func.get_inline_value(_name);
    return _current_val;
  }

  /**
   * Returns a new copy of the accessor.
   */
  virtual UniquePtr<ParameterAccessor<T> > clone() const {
    return UniquePtr<ParameterAccessor<T> >
      (new ParsedFunctionParameter<T>(_func, _name));
  }

private:
  ParsedFunction<T> & _func;
  std::string _name;

  // We need to return a reference from get().  That's a pointless
  // pessimization for libMesh::Number but it might become worthwhile
  // later when we handle field parameters.
  mutable libMesh::Number _current_val;
};

} // namespace libMesh

#endif // LIBMESH_PARSED_FUNCTION_PARAMETER_H
