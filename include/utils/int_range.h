// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h" // cast_int

// libMesh includes
#include "numeric_vector.h"

// C++ includes
#include <vector>

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
 * We perform a cast_int operation (no-op in opt mode, test+assert in
 * debug) at construction time to make sure that the given range
 * bounds are representable by the given range type.
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

  template <typename U, typename V>
  IntRange(U begin, V end) :
    _begin(cast_int<T>(begin)),
    _end(cast_int<T>(end))
  {}

  iterator begin() const { return _begin; }

  iterator end () const { return _end; }

private:
  iterator _begin, _end;
};



/**
 * Helper function that returns an IntRange<std::size_t> representing
 * all the indices of the passed-in vector.
 */
template <typename T>
IntRange<std::size_t> index_range(const std::vector<T> & vec)
{
  return IntRange<std::size_t>(0, vec.size());
}



/**
 * Same thing but for NumericVector. Returns a range (first_local_index, last_local_index).
 */
template <typename T>
IntRange<numeric_index_type> index_range(const NumericVector<T> & vec)
{
  return {vec.first_local_index(), vec.last_local_index()};
}

} // namespace libMesh

#endif // LIBMESH_INT_RANGE_H
