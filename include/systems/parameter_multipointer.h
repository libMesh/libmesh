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



#ifndef LIBMESH_PARAMETER_MULTIPOINTER_H
#define LIBMESH_PARAMETER_MULTIPOINTER_H


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
 * This is a slightly flexible ParameterAccessor subclass: it stores
 * all user-provided pointers to copies of the parameter, and modifies
 * the value at each location in memory.
 *
 * \author Roy Stogner
 * \date 2015
 * \brief Stores multiple user-provided pointers.
 */
template <typename T=Number>
class ParameterMultiPointer : public ParameterAccessor<T>
{
public:
  /**
   * Constructor: no parameters attached yet
   */
  ParameterMultiPointer() {}

  /**
   * Constructor: take the first raw pointer to the parameter
   */
  ParameterMultiPointer(T * param_ptr) : _ptrs(1, param_ptr) {}

  /**
   * A simple reseater won't work with a multipointer
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
    libmesh_assert(!_ptrs.empty());
#ifndef NDEBUG
    // Compare other values to the last one we'll change
    const T & val = *_ptrs.back();
#endif
    for (std::size_t i=0; i != _ptrs.size(); ++i)
      {
        // If you're already using inconsistent parameters we can't
        // help you.
        libmesh_assert_equal_to(*_ptrs[i], val);
        *_ptrs[i] = new_value;
      }
  }

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T & get () const libmesh_override
  {
    libmesh_assert(!_ptrs.empty());
    T & val = *_ptrs[0];
#ifndef NDEBUG
    // If you're already using inconsistent parameters we can't help
    // you.
    for (std::size_t i=1; i < _ptrs.size(); ++i)
      libmesh_assert_equal_to(*_ptrs[i], val);
#endif
    return val;
  }

  /**
   * Returns a new copy of the accessor.
   */
  virtual UniquePtr<ParameterAccessor<T> > clone() const libmesh_override
  {
    ParameterMultiPointer * pmp = new ParameterMultiPointer<T>();
    pmp->_ptrs = _ptrs;

    return UniquePtr<ParameterAccessor<T> >(pmp);
  }

  void push_back (T * new_ptr) { _ptrs.push_back(new_ptr); }

  /**
   * Returns the number of data associated with this parameter.
   * Useful for testing if the multipointer is empty/invalid.
   */
  std::size_t size() const { return _ptrs.size(); }

private:
  std::vector<T *> _ptrs;
};

} // namespace libMesh

#endif // LIBMESH_PARAMETER_MULTIPOINTER_H
