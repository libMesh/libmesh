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



#ifndef LIBMESH_INT_RANGE_H
#define LIBMESH_INT_RANGE_H

// C++ Includes   -----------------------------------
#include <map>

namespace libMesh
{

/**
 * The \p IntRange templated class is intended to make it easy to
 * loop over integers which are indices of a container.
 *
 * In cases where such a range is defined by the result of a virtual
 * function call, this allows range-based for loops to be easily
 * written which make only a single such call, rather than a new call
 * for each iteration.
 *
 * \author  Roy H. Stogner
 */

template <typename T>
class IntRange
{
public:
  class iterator {
  public:
    iterator (T i) : _i(i) {}

    T operator* () const { return _i; }

    const iterator & operator++ () {
      ++_i;
      return *this;
    }

    iterator operator++ (int) {
      iterator returnval(*this);
      ++_i;
      return returnval;
    }

    bool operator== (const iterator & j) const {
      return ( _i == j._i );
    }

    bool operator!= (const iterator & j) const {
      return !(*this == j);
    }

  private:
    T _i;
  };

  IntRange(T begin, T end) : _begin(begin), _end(end) {}

  iterator begin() const { return _begin; }

  iterator end () const { return _end; }

private:
  iterator _begin, _end;
};

} // namespace libMesh

#endif // LIBMESH_INT_RANGE_H
