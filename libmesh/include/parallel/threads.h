// $Id: parallel.h 2631 2008-02-04 16:46:11Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


#ifndef __threads_h__
#define __threads_h__

// System includes

// Local includes
#include "libmesh_config.h"

// Threading building blocks includes
#ifdef HAVE_TBB_API
#  include "tbb/tbb_stddef.h"
#  include "tbb/blocked_range.h"
#  include "tbb/parallel_for.h"
#  include "tbb/parallel_reduce.h"
#  include "tbb/task_scheduler_init.h"
#  include "tbb/spin_mutex.h"
#  include "tbb/atomic.h"
#endif

/**
 * The Threads namespace is for wrapper functions
 * for common general multithreading algorithms and tasks.
 */
namespace Threads
{
#ifdef HAVE_TBB_API
  //-------------------------------------------------------------------
  /**
   * Scheduler to manage threads.
   */
  typedef tbb::task_scheduler_init task_scheduler_init;

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
  { tbb::parallel_for (range, body); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified 
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_for (const Range &range, const Body &body, const Partitioner &partitioner)
  { tbb::parallel_for (range, body, partitioner); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified 
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_reduce (const Range &range, Body &body)
  { tbb::parallel_reduce (range, body); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified 
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_reduce (const Range &range, Body &body, const Partitioner &partitioner)
  { tbb::parallel_reduce (range, body, partitioner); }

  //-------------------------------------------------------------------
  /**
   * Spin mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.
   */
  typedef tbb::spin_mutex spin_mutex;

  //-------------------------------------------------------------------
  /**
   * Defines atomic operations which can only be executed on a 
   * single thread at a time.  This is used in reference counting,
   * for example, to allow count++/count-- to work.
   */
  template <typename T>
  class atomic : public tbb::atomic<T> {};

  

#else //HAVE_TBB_API

  //-------------------------------------------------------------------
  /**
   * Scheduler to manage threads.
   */
  class task_scheduler_init
  {
  public:
    static const int automatic = -1;
    task_scheduler_init (int = automatic) {};
    void initialize (int = automatic) {};
    void terminate () {};    
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
  { body(range); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided function object in parallel on the specified 
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_for (const Range &range, const Body &body, const Partitioner &)
  { body(range); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified 
   * range.
   */
  template <typename Range, typename Body>
  inline
  void parallel_reduce (const Range &range, Body &body)
  { body(range); }

  //-------------------------------------------------------------------
  /**
   * Exectue the provided reduction operation in parallel on the specified 
   * range with the specified partitioner.
   */
  template <typename Range, typename Body, typename Partitioner>
  inline
  void parallel_reduce (const Range &range, Body &body, const Partitioner &)
  { body(range); }

  //-------------------------------------------------------------------
  /**
   * Spin mutex.  Implements mutual exclusion by busy-waiting in user
   * space for the lock to be acquired.
   */
  class spin_mutex 
  {
  public:
    spin_mutex() {}

    class scoped_lock 
    {
    public:
      scoped_lock () {}
      scoped_lock ( spin_mutex&  ) {} 
      void acquire ( spin_mutex& ) {}
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

  

#endif // #ifdef HAVE_TBB_API

  /**
   * A spin mutex object which 
   */
  extern spin_mutex spin_mtx;



} // namespace Threads

#endif // #define __threads_h__
