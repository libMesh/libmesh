// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_THREADS_H
#define LIBMESH_THREADS_H

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"  // for libmesh_assert


// Compile-time check: TBB and pthreads are now mutually exclusive.
#if defined(LIBMESH_HAVE_TBB_API) && defined(LIBMESH_HAVE_PTHREAD)
MULTIPLE THREADING MODELS CANNOT BE SIMULTANEOUSLY ACTIVE
#endif

namespace libMesh
{

/**
 * The Threads namespace is for wrapper functions
 * for common general multithreading algorithms and tasks.
 */
namespace Threads
{

/**
 * A boolean which is true iff we are in a Threads:: function
 * It may be useful to assert(!Threads::in_threads) in any code
 * which is known to not be thread-safe.
 */
extern bool in_threads;

/**
 * We use a class to turn Threads::in_threads on and off, to be
 * exception-safe.
 */
class BoolAcquire
{
public:
  explicit
  BoolAcquire(bool & b) : _b(b) { libmesh_assert(!_b); _b = true; }

  ~BoolAcquire() { libmesh_exceptionless_assert(_b); _b = false; }
private:
  bool & _b;
};


/**
 * Simple compatibility class for std::thread 'concurrent' execution.
 * Not at all concurrent, but provides a compatible interface.
 */
class NonConcurrentThread
{
public:
  /**
   * Constructor.  Takes a callable function object and executes it.
   * Our wrapper class actually blocks execution until the thread
   * is complete.
   */
  template <typename Callable>
  NonConcurrentThread (Callable f) { f(); }

  /**
   * Join is a no-op, since the constructor blocked until completion.
   */
  void join() {}

  /**
   * Always joinable.
   */
  bool joinable() const { return true; }
};

} // namespace Threads

} // namespace libMesh



// Include thread-model specific algorithms and objects.  These
// headers include headers of their own and handle their own
// namespacing.
#ifdef LIBMESH_HAVE_TBB_API
# include "libmesh/threads_tbb.h"
#elif LIBMESH_HAVE_PTHREAD
# include "libmesh/threads_pthread.h"
#else
# include "libmesh/threads_none.h"
#endif



namespace libMesh
{

namespace Threads
{

/**
 * Blocked range which can be subdivided and executed in parallel.
 */
template <typename T>
class BlockedRange
{
public:
  /**
   * Allows an \p StoredRange to behave like an STL container.
   */
  typedef T const_iterator;

  /**
   * Constructor. Optionally takes the \p grainsize parameter, which is the
   * smallest chunk the range may be broken into for parallel
   * execution.
   */
  explicit BlockedRange (const unsigned int new_grainsize = 1000) :
    _grainsize(new_grainsize)
  {}

  /**
   * Constructor.  Takes the beginning and end of the range.
   * Optionally takes the \p grainsize parameter, which is the
   * smallest chunk the range may be broken into for parallel
   * execution.
   */
  BlockedRange (const const_iterator first,
                const const_iterator last,
                const unsigned int new_grainsize = 1000) :
    _grainsize(new_grainsize)
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
  BlockedRange (const BlockedRange<T> & r):
    _end(r._end),
    _begin(r._begin),
    _grainsize(r._grainsize)
  {}

  /**
   * Splits the range \p r.  The first half
   * of the range is left in place, the second
   * half of the range is placed in *this.
   */
  BlockedRange (BlockedRange<T> & r, Threads::split ) :
    _end(r._end),
    _begin(r._begin),
    _grainsize(r._grainsize)
  {
    const_iterator
      beginning = r._begin,
      ending    = r._end,
      middle    = beginning + (ending - beginning)/2u;

    r._end = _begin = middle;
  }

  /**
   * Resets the \p StoredRange to contain [first,last).
   */
  void reset (const const_iterator first,
              const const_iterator last)
  {
    _begin = first;
    _end   = last;
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
   * The grain size for the range.  The range will be subdivided into
   * subranges not to exceed the grain size.
   */
  unsigned int grainsize () const {return _grainsize;}

  /**
   * Set the grain size.
   */
  void grainsize (const unsigned int & gs) {_grainsize = gs;}

  /**
   * \return the size of the range.
   */
  int size () const { return (_end -_begin); }

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
  bool is_divisible() const { return ((_begin + this->grainsize()) < _end); }

private:

  const_iterator _end;
  const_iterator _begin;
  unsigned int _grainsize;
};



/**
 * A convenient spin mutex object which can be used for obtaining locks.
 */
extern spin_mutex spin_mtx;

/**
 * A convenient recursive mutex object which can be used for obtaining locks.
 */
extern recursive_mutex recursive_mtx;

} // namespace Threads

} // namespace libMesh

#endif // LIBMESH_THREADS_H
