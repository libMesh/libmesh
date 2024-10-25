// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes
#include <vector>

namespace libMesh
{

// Forward declarations
template <typename T> class DenseSubVector;
template <typename T> class DenseVector;
template <typename T> class NumericVector;

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
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = T;
    using pointer = T;
    using reference = T&;

    iterator (T i) : _i(i) {}

    T operator* () const { return _i; }

    const iterator & operator++ ()
    {
      ++_i;
      return *this;
    }

    iterator operator++ (int)
    {
      iterator returnval(*this);
      ++_i;
      return returnval;
    }

    bool operator== (const iterator & j) const
    {
      return ( _i == j._i );
    }

    bool operator!= (const iterator & j) const
    {
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
 * all the indices of the passed-in vector-like object (i.e. any type that
 * has a size() member).
 */
template <typename T>
auto index_range(const T & sizable)
{
  return IntRange<decltype(sizable.size())>(0, sizable.size());
}


/**
 * Same thing but for NumericVector. Returns a range (first_local_index, last_local_index).
 */
template <typename T>
IntRange<numeric_index_type> index_range(const NumericVector<T> & vec)
{
  return {vec.first_local_index(), vec.last_local_index()};
}



/**
 * The 2-parameter make_range() helper function returns an IntRange<T>
 * when both input parameters are of type T. This saves a bit of
 * typing over calling the IntRange<T> constructor directly.
 */
template <typename T>
IntRange<T> make_range(T beg, T end)
{
  return {beg, end};
}



/**
 * The 1-parameter version of make_range() saves even more typing in
 * the common case of a 0 starting point. Example usage:
 *
 * for (auto i : make_range(10))
 *
 * will loop from 0 to 9. In more realistic cases such as:
 *
 * for (auto i : make_range(foo()))
 *
 * this construction guarantees that the function foo() is called
 * exactly once rather than once per loop iteration.
 */
template <typename T>
IntRange<T> make_range(T end)
{
  return {T(0), end};
}

} // namespace libMesh

#endif // LIBMESH_INT_RANGE_H
