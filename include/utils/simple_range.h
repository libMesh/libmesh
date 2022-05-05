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



#ifndef LIBMESH_SIMPLERANGE_H
#define LIBMESH_SIMPLERANGE_H

#include <utility>

namespace libMesh
{

/**
 * The \p SimpleRange templated class is intended to make it easy to
 * construct ranges from pairs of iterators.
 *
 * \author  Roy H. Stogner
 */

template <typename IndexType>
class SimpleRange
{
public:
  SimpleRange(IndexType begin, IndexType end) : _begin(begin), _end(end) {}

  IndexType begin() const { return _begin; }

  IndexType end () const { return _end; }

private:
  IndexType _begin, _end;
};



/**
 * Helper function that allows us to treat a homogenous pair as a
 * range. Useful for writing range-based for loops over the pair
 * returned by std::equal_range() and std::map::equal_range().
 */
template<typename IndexType>
SimpleRange<IndexType> as_range(const std::pair<IndexType, IndexType> & p)
{
  return {p.first, p.second};
}



/**
 * As above, but can be used in cases where a std::pair is not
 * otherwise involved.
 */
template<typename IndexType>
SimpleRange<IndexType> as_range(const IndexType & first,
                        const IndexType & second)
{
  return {first, second};
}

} // namespace libMesh

#endif // LIBMESH_SIMPLERANGE_H
