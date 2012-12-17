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


#ifndef LIBMESH_THREADS_H
#define LIBMESH_THREADS_H

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"  // for libmesh_assert

// Threading building blocks includes
#ifdef LIBMESH_HAVE_TBB_API
#  include "libmesh/libmesh_logging.h" // only mess with the perflog if we are really multithreaded
#  include "tbb/tbb_stddef.h"
#  include "tbb/blocked_range.h"
#  include "tbb/parallel_for.h"
#  include "tbb/parallel_reduce.h"
#  include "tbb/task_scheduler_init.h"
#  include "tbb/partitioner.h"
#  include "tbb/spin_mutex.h"
#  include "tbb/recursive_mutex.h"
#  include "tbb/atomic.h"
#endif

// C++ includes
#ifdef LIBMESH_HAVE_STD_THREAD
#  include <thread>
#elif LIBMESH_HAVE_TBB_CXX_THREAD
#  include "tbb/tbb_thread.h"
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
  class BoolAcquire {
  public:
    explicit
    BoolAcquire(bool& b) : _b(b) { libmesh_assert(!_b); _b = true; }

    ~BoolAcquire() { libmesh_assert(_b); _b = false; }
  private:
    bool& _b;
  };


#ifdef LIBMESH_HAVE_STD_THREAD
  //--------------------------------------------------------------------
  /**
   * Use std::thread when available.
   */
  typedef std::thread Thread;

#elif LIBMESH_HAVE_TBB_CXX_THREAD
  //--------------------------------------------------------------------
  /**
   * Fall back to tbb::tbb_thread when available.
   */
  typedef tbb::tbb_thread Thread;
  
#else
  //--------------------------------------------------------------------
  /**
   * Simple compatibility class for std::thread 'concurrent' execution.
   * Not at all concurrent, but provides a compatible interface.
   */
  class Thread
  {
  public:
    /**
     * Constructor.  Takes a callable function object and executes it.
     * Our wrapper class actually blocks execution until the thread
     * is complete.
     */
    template <typename Callable>
    Thread (Callable f) { f(); }

    /**
     * Join is a no-op, since the constructor blocked until completion.
     */
    void join() {}

    /**
     * Always joinable.
     */
    bool joinable() const { return true; }
  };
  
#endif



#ifdef LIBMESH_HAVE_TBB_API
  //-------------------------------------------------------------------
  /**
   * Scheduler to manage threads.
   */
  typedef tbb::task_scheduler_init task_scheduler_init;



  //-------------------------------------------------------------------
  /**
   * Dummy "splitting object" used to distinguish splitting constructors
   * from copy constructors.
   */
  typedef tbb::split split;



  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_for (const Range &range, const Body &body)
  {
    BoolAcquire b(in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    const bool logging_was_enabled = libMesh::perflog.logging_enabled();

    if (libMesh::n_threads() > 1)
      libMesh::perflog.disable_logging();
#endif

    if (libMesh::n_threads() > 1)
      tbb::parallel_for (range, body, tbb::auto_partitioner());

    else
      body(range);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    if (libMesh::n_threads() > 1 && logging_was_enabled)
      libMesh::perflog.enable_logging();
#endif
  }



  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_for (const Range &range, const Body &body, const Partitioner &partitioner)
  {
    BoolAcquire b(in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    const bool logging_was_enabled = libMesh::perflog.logging_enabled();

    if (libMesh::n_threads() > 1)
      libMesh::perflog.disable_logging();
#endif

    if (libMesh::n_threads() > 1)
      tbb::parallel_for (range, body, partitioner);

    else
      body(range);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    if (libMesh::n_threads() > 1 && logging_was_enabled)
      libMesh::perflog.enable_logging();
#endif
  }



  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_reduce (const Range &range, Body &body)
  {
    BoolAcquire b(in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    const bool logging_was_enabled = libMesh::perflog.logging_enabled();

    if (libMesh::n_threads() > 1)
      libMesh::perflog.disable_logging();
#endif

    if (libMesh::n_threads() > 1)
      tbb::parallel_reduce (range, body, tbb::auto_partitioner());

    else
      body(range);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    if (libMesh::n_threads() > 1 && logging_was_enabled)
      libMesh::perflog.enable_logging();
#endif
  }



  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_reduce (const Range &range, Body &body, const Partitioner &partitioner)
  {
    BoolAcquire b(in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    const bool logging_was_enabled = libMesh::perflog.logging_enabled();

    if (libMesh::n_threads() > 1)
      libMesh::perflog.disable_logging();
#endif

       if (libMesh::n_threads() > 1)
	 tbb::parallel_reduce (range, body, partitioner);

       else
	 body(range);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
    if (libMesh::n_threads() > 1 && logging_was_enabled)
      libMesh::perflog.enable_logging();
#endif
  }



  //-------------------------------------------------------------------
  /**
   * Spin mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.
   */
  typedef tbb::spin_mutex spin_mutex;

  //-------------------------------------------------------------------
  /**
   * Recursive mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.  The same thread can aquire the
   * same lock multiple times
   */
  typedef tbb::recursive_mutex recursive_mutex;

  //-------------------------------------------------------------------
  /**
   * Defines atomic operations which can only be executed on a
   * single thread at a time.  This is used in reference counting,
   * for example, to allow count++/count-- to work.
   */
  template <typename T>
  class atomic : public tbb::atomic<T> {};



#else //LIBMESH_HAVE_TBB_API

  //-------------------------------------------------------------------
  /**
   * Scheduler to manage threads.
   */
  class task_scheduler_init
  {
  public:
    static const int automatic = -1;
    explicit task_scheduler_init (int = automatic) {}
    void initialize (int = automatic) {}
    void terminate () {}
  };

  //-------------------------------------------------------------------
  /**
   * Dummy "splitting object" used to distinguish splitting constructors
   * from copy constructors.
   */
  class split {};

  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_for (const Range &range, const Body &body)
  {
    BoolAcquire b(in_threads);
    body(range);
  }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_for (const Range &range, const Body &body, const Partitioner &)
  {
    BoolAcquire b(in_threads);
    body(range);
  }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_reduce (const Range &range, Body &body)
  {
    BoolAcquire b(in_threads);
    body(range);
  }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_reduce (const Range &range, Body &body, const Partitioner &)
  {
    BoolAcquire b(in_threads);
    body(range);
  }

  //-------------------------------------------------------------------
  /**
   * Spin mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.
   */
  class spin_mutex
  {
  public:
    spin_mutex() {}
    void lock () {}
    void unlock () {}

    class scoped_lock
    {
    public:
      scoped_lock () {}
      explicit scoped_lock ( spin_mutex&  ) {}
      void acquire ( spin_mutex& ) {}
      void release () {}
    };
  };

  //-------------------------------------------------------------------
  /**
   * Recursive mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.
   */
  class recursive_mutex
  {
  public:
    recursive_mutex() {}

    class scoped_lock
    {
    public:
      scoped_lock () {}
      explicit scoped_lock ( recursive_mutex&  ) {}
      void acquire ( recursive_mutex& ) {}
      void release () {}
    };
  };

  //-------------------------------------------------------------------
  /**
   * Defines atomic operations which can only be executed on a
   * single thread at a time.
   */
  template <typename T>
  class atomic
  {
  public:
    atomic () : _val(0) {}
    operator T& () { return _val; }
  private:
    T _val;
  };


#endif // #ifdef LIBMESH_HAVE_TBB_API



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
    explicit BlockedRange (const unsigned int grainsize = 1000) :
      _grainsize(grainsize)
    {}

    /**
     * Constructor.  Takes the beginning and end of the range.
     * Optionally takes the \p grainsize parameter, which is the
     * smallest chunk the range may be broken into for parallel
     * execution.
     */
    BlockedRange (const const_iterator first,
		  const const_iterator last,
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
    BlockedRange (const BlockedRange<T> &r):
      _end(r._end),
      _begin(r._begin),
      _grainsize(r._grainsize)
    {}

    /**
     * Splits the range \p r.  The first half
     * of the range is left in place, the second
     * half of the range is placed in *this.
     */
    BlockedRange (BlockedRange<T> &r, Threads::split ) :
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
    void grainsize (const unsigned int &gs) {_grainsize = gs;}

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
   * A spin mutex object which
   */
  extern spin_mutex spin_mtx;

  /**
   * A recursive mutex object which
   */
  extern recursive_mutex recursive_mtx;

} // namespace Threads

} // namespace libMesh

#endif // LIBMESH_THREADS_H
