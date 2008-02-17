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
#  include "libmesh_logging.h" // only mess with the perflog if we are really multithreaded
#  include "tbb/tbb_stddef.h"
#  include "tbb/blocked_range.h"
#  include "tbb/parallel_for.h"
#  include "tbb/parallel_reduce.h"
#  include "tbb/task_scheduler_init.h"
#  include "tbb/partitioner.h"
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
#ifdef ENABLE_PERFORMANCE_LOGGING
    libMesh::perflog.disable_logging();
#endif   

    tbb::parallel_for (range, body, tbb::auto_partitioner()); 

#ifdef ENABLE_PERFORMANCE_LOGGING
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
#ifdef ENABLE_PERFORMANCE_LOGGING
    libMesh::perflog.disable_logging();
#endif   

    tbb::parallel_for (range, body, partitioner); 

#ifdef ENABLE_PERFORMANCE_LOGGING
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
#ifdef ENABLE_PERFORMANCE_LOGGING
    libMesh::perflog.disable_logging();
#endif   

    tbb::parallel_reduce (range, body, tbb::auto_partitioner()); 

#ifdef ENABLE_PERFORMANCE_LOGGING
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
#ifdef ENABLE_PERFORMANCE_LOGGING
    libMesh::perflog.disable_logging();
#endif   

    tbb::parallel_reduce (range, body); 

#ifdef ENABLE_PERFORMANCE_LOGGING
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



//   /**
//    * Blocked range which can be subdivided and executed in parallel.
//    */
//   template <typename T>
//   class BlockedRange 
//   {
//   public:
//     /**
//      * Allows an \p StoredRange to behave like an STL container.
//      */
//     typedef T const_iterator;

//     /**
//      * Constructor. Optionally takes the \p grainsize parameter, which is the
//      * smallest chunk the range may be broken into for parallel
//      * execution.
//      */
//     BlockedRange (const unsigned int grainsize = 1000) :
//       _grainsize(grainsize)
//     {}

//     /**
//      * Constructor.  Takes the beginning and end of the range.  
//      * Optionally takes the \p grainsize parameter, which is the
//      * smallest chunk the range may be broken into for parallel
//      * execution.
//      */
//     BlockedRange (const const_iterator first,
// 		  const const_iterator last,
// 		  const unsigned int grainsize = 1000) :
//       _grainsize(grainsize)
//     {
//       this->reset(first, last);
//     }

//     /**
//      * Copy constructor.  The \p StoredRange can be copied into
//      * subranges for parallel execution.  In this way the 
//      * initial \p StoredRange can be thought of as the root of
//      * a binary tree.  The root element is the only element
//      * which interacts with the user.  It takes a specified
//      * range of objects and packs it into a contiguous vector
//      * which can be split efficiently. However, there is no need
//      * for the child ranges to contain this vector, so long as
//      * the parent outlives the children.  So we implement
//      * the copy constructor to specifically omit the \p _objs
//      * vector.
//      */
//     BlockedRange (const BlockedRange<T> &r):
//       _end(r._end), 
//       _begin(r._begin),
//       _grainsize(r._grainsize)
//     {}
    
//     /**
//      * Splits the range \p r.  The first half
//      * of the range is left in place, the second
//      * half of the range is placed in *this.
//      */
//     BlockedRange (BlockedRange<T> &r, Threads::split ) : 
//       _end(r._end),
//       _begin(r._begin),
//       _grainsize(r._grainsize)
//     {
//       const_iterator
// 	beginning = r._begin,
// 	ending    = r._end,
// 	middle    = beginning + (ending - beginning)/2u;
    
//       r._end = _begin = middle;
//     }
    
//     /**
//      * Resets the \p StoredRange to contain [first,last).
//      */ 
//     void reset (const const_iterator first,
// 		const const_iterator last)
//     {
//       _begin = first;
//       _end   = last;
//     }
    
//     /**
//      * Beginning of the range.
//      */
//     const_iterator begin () const { return _begin; }  
    
//     /**
//      * End of the range.
//      */
//     const_iterator end () const { return _end; }
    
//     /**
//      * The grain size for the range.  The range will be subdivided into 
//      * subranges not to exceed the grain size.
//      */
//     unsigned int grainsize () const {return _grainsize;}
    
//     /**
//      * Set the grain size.
//      */
//     void grainsize (const unsigned int &gs) {_grainsize = gs;}
    
//     /**
//      * \return the size of the range.
//      */
//     int size () const { return (_end -_begin); }

//     //------------------------------------------------------------------------
//     // Methods that implement Range concept
//     //------------------------------------------------------------------------
    
//     /**
//      * Returns true if the range is empty.
//      */  
//     bool empty() const { return (_begin == _end); }
    
//     /**
//      * Returns true if the range can be subdivided.
//      */
//     bool is_divisible() const { return static_cast<int>(this->grainsize()) < (_end - _begin); }

//   private:
    
//     const_iterator _end;
//     const_iterator _begin;
//     unsigned int _grainsize;
//   };


  
  /**
   * A spin mutex object which 
   */
  extern spin_mutex spin_mtx;



} // namespace Threads

#endif // #define __threads_h__
