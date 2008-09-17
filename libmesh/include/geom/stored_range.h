// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __stored_range_h__
#define __stored_range_h__

// C++ includes
#include <vector>

// Local includes
#include "threads.h"


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
  StoredRange (const unsigned int grainsize = 1000) :
  _grainsize(grainsize)
  {}

  /**
   * Constructor.  Takes the beginning and end of the range.  
   * Optionally takes the \p grainsize parameter, which is the
   * smallest chunk the range may be broken into for parallel
   * execution.
   */
  StoredRange (const iterator_type &first,
	       const iterator_type &last,
	       const unsigned int grainsize = 1000) :
  _grainsize(grainsize)
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
    _grainsize(er._grainsize)
  {
    // specifically, do *not* copy the vector
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
    _grainsize(r._grainsize)
  {
    const_iterator
      beginning = r._begin,
      ending    = r._end,
      middle    = beginning + std::distance(beginning, ending)/2u;
    
    r._end = _begin = middle;

    unsigned int
      first = r._first,
      last  = r._last,
      half  = first + (last-first)/2u;

    r._last = _first = half;
  }
    
  /**
   * Resets the \p StoredRange to contain [first,last).
   */ 
  void reset (const iterator_type &first,
	      const iterator_type &last)
  {
    _objs.clear();

    for (iterator_type it=first; it!=last; ++it)
      _objs.push_back(*it);
    
    _begin = _objs.begin();
    _end   = _objs.end();

    _first = 0;
    _last  = _objs.size();
  }

  /**
   * Resets the range to the last specified range.  This method only exists
   * for efficiency -- it is more efficient to set the range to its previous
   * value without rebuilding the underlying vector.
   */
  void reset ()
  {
    _begin = _objs.begin();
    _end   = _objs.end();

    _first = 0;
    _last  = _objs.size();
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
  unsigned int first_idx () const { return _first; }

  /**
   * Index in the stored vector of the last object.
   */
  unsigned int last_idx () const { return _last; }

  /**
   * The grain size for the range.  The range will be subdivided into 
   * subranges not to exceed the grain size.
   */
  unsigned int grainsize () const {return _grainsize;}

  /**
   * Set the grain size.
   */
  void  grainsize (const unsigned int &gs) {_grainsize = gs;}

  /**
   * \return the size of the range.
   */
  unsigned int size () const { return std::distance(_begin, _end); }

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
  unsigned int _last;
  unsigned int _first;
  unsigned int _grainsize;
  std::vector<object_type> _objs;
};

#endif // end #ifndef __stored_range_h__
