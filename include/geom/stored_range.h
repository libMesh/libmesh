// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_STORED_RANGE_H
#define LIBMESH_STORED_RANGE_H

// Local includes
#include "libmesh/threads.h"

// C++ includes
#include <vector>

namespace libMesh
{


/**
 * The \p StoredRange class defined a contiguous, divisible set of objects
 * This class is used primarily as the argument to function objects.  The
 * range can then be subdivided into a number of "tasks" which can be executed
 * in parallel.  This concept is central to the Threading Building Blocks
 * template library which can optionally be used by \p libMesh to implement
 * shared-memory parallelism.
 *
 * The implementation takes a user-provided object range and packs it into
 * a contiguous vector which can then be subdivided efficiently.  A first-cut
 * implementation using raw element iterators incurred simply too much overhead
 * by using the predicated iterators, specifically operations such as advancing
 * such iterators has a cost proportional to the amount the iterator is advanced.
 * Hence in this implementation the user-provided range is packed into a vector.
 *
 * \author Benjamin S. Kirk, 2008.
 */
template <typename iterator_type, typename object_type>
class StoredRange
{
public:
  /**
   * Allows an \p StoredRange to behave like an STL container.
   */
  typedef typename std::vector<object_type>::const_iterator const_iterator;

  /**
   * Constructor. Optionally takes the \p grainsize parameter, which is the
   * smallest chunk the range may be broken into for parallel
   * execution.
   */
  StoredRange (const unsigned int new_grainsize = 1000) :
    _end(),
    _begin(),
    _last(),
    _first(),
    _grainsize(new_grainsize),
    _objs()
  {}

  /**
   * Constructor.  Takes the beginning and end of the range.
   * Optionally takes the \p grainsize parameter, which is the
   * smallest chunk the range may be broken into for parallel
   * execution.
   */
  StoredRange (const iterator_type &first,
	       const iterator_type &last,
	       const unsigned int new_grainsize = 1000) :
    _end(),
    _begin(),
    _last(),
    _first(),
    _grainsize(new_grainsize),
    _objs()
  {
    this->reset(first, last);
  }

  /**
   * Copy constructor.  The \p StoredRange can be copied into
   * subranges for parallel execution.  In this way the
   * initial \p StoredRange can be thought of as the root of
   * a binary tree.  The root element is the only element
   * which interacts with the user.  It takes a specified
   * range of objects and packs it into a contiguous vector
   * which can be split efficiently. However, there is no need
   * for the child ranges to contain this vector, so long as
   * the parent outlives the children.  So we implement
   * the copy constructor to specifically omit the \p _objs
   * vector.
   */
  StoredRange (const StoredRange<iterator_type,object_type> &er):
    _end(er._end),
    _begin(er._begin),
    _last(er._last),
    _first(er._first),
    _grainsize(er._grainsize),
    _objs()
  {
    // specifically, do *not* copy the vector
  }

  /**
   * NOTE: When using pthreads this constructor is MANDATORY!!!
   *
   * Copy constructor.  The \p StoredRange can be copied into
   * subranges for parallel execution.  In this way the
   * initial \p StoredRange can be thought of as the root of
   * a binary tree.  The root element is the only element
   * which interacts with the user.  It takes a specified
   * range of objects and packs it into a contiguous vector
   * which can be split efficiently. However, there is no need
   * for the child ranges to contain this vector, so long as
   * the parent outlives the children.  So we implement
   * the copy constructor to specifically omit the \p _objs
   * vector. This version allows you to set the beginning and
   * ending of this new range to be different from that of the
   * one we're copying.
   */
  StoredRange (const StoredRange<iterator_type,object_type> &er,
               const const_iterator &begin_range,
               const const_iterator &end_range):
    _end(end_range),
    _begin(begin_range),
    _last(0), // Initialize these in a moment
    _first(0),
    _grainsize(er._grainsize),
    _objs()
  {
    // specifically, do *not* copy the vector

    _first = std::distance(er._begin, _begin);
    _last = _first + std::distance(_begin, _end);
  }

  /**
   * Splits the range \p r.  The first half
   * of the range is left in place, the second
   * half of the range is placed in *this.
   */
  StoredRange (StoredRange<iterator_type,object_type> &r, Threads::split ) :
    _end(r._end),
    _begin(r._begin),
    _last(r._last),
    _first(r._first),
    _grainsize(r._grainsize),
    _objs()
  {
    const_iterator
      beginning = r._begin,
      ending    = r._end,
      middle    = beginning + std::distance(beginning, ending)/2u;

    r._end = _begin = middle;

    std::size_t
      first = r._first,
      last  = r._last,
      half  = first + (last-first)/2u;

    r._last = _first = half;
  }

  /**
   * Resets the \p StoredRange to contain [first,last).  Returns
   * a reference to itself for convenience, so functions
   * expecting a StoredRange<> can be passed e.g. foo.reset(begin,end).
   */
  StoredRange<iterator_type, object_type> &
  reset (const iterator_type &first,
	 const iterator_type &last)
  {
    _objs.clear();

    for (iterator_type it=first; it!=last; ++it)
      _objs.push_back(*it);

    _begin = _objs.begin();
    _end   = _objs.end();

    _first = 0;
    _last  = _objs.size();

    return *this;
  }

  /**
   * Resets the range to the last specified range.  This method only exists
   * for efficiency -- it is more efficient to set the range to its previous
   * value without rebuilding the underlying vector.  Returns
   * a reference to itself for convenience, so functions
   * expecting a StoredRange<> can be passed e.g. foo.reset().
   */
  StoredRange<iterator_type, object_type> & reset ()
  {
    _begin = _objs.begin();
    _end   = _objs.end();

    _first = 0;
    _last  = _objs.size();

    return *this;
  }

  /**
   * Beginning of the range.
   */
  const_iterator begin () const { return _begin; }

  /**
   * End of the range.
   */
  const_iterator end () const { return _end; }

  /**
   * Index in the stored vector of the first object.
   */
  std::size_t first_idx () const { return _first; }

  /**
   * Index in the stored vector of the last object.
   */
  std::size_t last_idx () const { return _last; }

  /**
   * The grain size for the range.  The range will be subdivided into
   * subranges not to exceed the grain size.
   */
  std::size_t grainsize () const {return _grainsize;}

  /**
   * Set the grain size.
   */
  void grainsize (const unsigned int &gs) {_grainsize = gs;}

  /**
   * \return the size of the range.
   */
  std::size_t size () const { return std::distance(_begin, _end); }

  //------------------------------------------------------------------------
  // Methods that implement Range concept
  //------------------------------------------------------------------------

  /**
   * Returns true if the range is empty.
   */
  bool empty() const { return (_begin == _end); }

  /**
   * Returns true if the range can be subdivided.
   */
  bool is_divisible() const { return this->grainsize() < static_cast<unsigned int>(std::distance(_begin, _end)); }

private:

  const_iterator _end;
  const_iterator _begin;
  std::size_t _last;
  std::size_t _first;
  std::size_t _grainsize;
  std::vector<object_type> _objs;
};

} // namespace libMesh

#endif // LIBMESH_STORED_RANGE_H
