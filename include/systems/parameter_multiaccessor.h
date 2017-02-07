// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARAMETER_MULTIACCESSOR_H
#define LIBMESH_PARAMETER_MULTIACCESSOR_H


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_accessor.h"

// C++ Includes
#include <vector>

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
 * \brief Stores a user-provided pointer to a parameter.
 */
template <typename T=Number>
class ParameterMultiAccessor : public ParameterAccessor<T>
{
public:
  /**
   * Constructor: no parameters attached yet
   */
  ParameterMultiAccessor() {}

  /**
   * Constructor: take the first sub-accessor for the parameter
   */
  ParameterMultiAccessor(const ParameterAccessor<T> & param_accessor) :
    _accessors(1, param_accessor.clone()) {}

  /*
   * Destructor: delete our clones of sub-accessors
   */
  ~ParameterMultiAccessor() {
    for (std::size_t i=0; i != _accessors.size(); ++i)
      delete _accessors[i];
  }

  /**
   * A simple reseater won't work with a multi-accessor
   */
  virtual ParameterAccessor<T> &
  operator= (T * /* new_ptr */) libmesh_override
  {
    libmesh_error();
    return *this;
  }

  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const T & new_value) libmesh_override
  {
    libmesh_assert(!_accessors.empty());
#ifndef NDEBUG
    // Compare other values to the last one we'll change
    const T & val = _accessors.back()->get();
#endif
    for (std::size_t i=0; i != _accessors.size(); ++i)
      {
        // If you're already using inconsistent parameters we can't
        // help you.
        libmesh_assert_equal_to(_accessors[i]->get(), val);
        _accessors[i]->set(new_value);
      }
  }

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const libmesh_override
  {
    libmesh_assert(!_accessors.empty());
    const T & val = _accessors[0]->get();
#ifndef NDEBUG
    // If you're already using inconsistent parameters we can't help
    // you.
    for (std::size_t i=1; i < _accessors.size(); ++i)
      libmesh_assert_equal_to(_accessors[i]->get(), val);
#endif
    return val;
  }

  /**
   * Returns a new copy of the accessor.
   */
  virtual UniquePtr<ParameterAccessor<T> > clone() const libmesh_override
  {
    ParameterMultiAccessor * pmp = new ParameterMultiAccessor<T>();
    for (std::size_t i=0; i != _accessors.size(); ++i)
      pmp->_accessors.push_back(_accessors[i]->clone().release());

    return UniquePtr<ParameterAccessor<T> >(pmp);
  }


  void push_back (const ParameterAccessor<T> & new_accessor) {
    _accessors.push_back(new_accessor.clone().release());
  }

  /**
   * Returns the number of sub-accessors associated with this
   * parameter.  Useful for testing if the multi-accessor is
   * empty/invalid.
   */
  std::size_t size() const { return _accessors.size(); }

private:
  std::vector<ParameterAccessor<T> *> _accessors;
};

} // namespace libMesh

#endif // LIBMESH_PARAMETER_MULTIACCESSOR_H
