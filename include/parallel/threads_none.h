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


#ifndef LIBMESH_THREADS_NONE_H
#define LIBMESH_THREADS_NONE_H

// Do not try to #include this header directly, it is designed to be
// #included directly by threads.h
#if !defined(LIBMESH_HAVE_TBB_API) && !defined(LIBMESH_HAVE_PTHREAD)

// Thread-Local-Storage macros
#define LIBMESH_TLS_TYPE(type)  type
#define LIBMESH_TLS_REF(value)  (value)

namespace libMesh
{

namespace Threads
{

/**
 * Use the non-concurrent placeholder.
 */
typedef NonConcurrentThread Thread;

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



/**
 * Dummy "splitting object" used to distinguish splitting constructors
 * from copy constructors.
 */
class split {};



/**
 * Execute the provided function object in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_for (const Range & range, const Body & body)
{
  BoolAcquire b(in_threads);
  body(range);
}



/**
 * Execute the provided function object in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_for (const Range & range, const Body & body, const Partitioner &)
{
  BoolAcquire b(in_threads);
  body(range);
}



/**
 * Execute the provided reduction operation in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_reduce (const Range & range, Body & body)
{
  BoolAcquire b(in_threads);
  body(range);
}



/**
 * Execute the provided reduction operation in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_reduce (const Range & range, Body & body, const Partitioner &)
{
  BoolAcquire b(in_threads);
  body(range);
}



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
    explicit scoped_lock ( spin_mutex &  ) {}
    void acquire ( spin_mutex & ) {}
    void release () {}
  };
};



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
    explicit scoped_lock ( recursive_mutex &  ) {}
    void acquire ( recursive_mutex & ) {}
    void release () {}
  };
};



/**
 * Defines atomic operations which can only be executed on a
 * single thread at a time.
 */
template <typename T>
class atomic
{
public:
  atomic () : _val(0) {}
  operator T & () { return _val; }
private:
  T _val;
};

} // namespace Threads

} // namespace libMesh

#endif // !defined(LIBMESH_HAVE_TBB_API) && !defined(LIBMESH_HAVE_PTHREAD)

#endif // LIBMESH_THREADS_NONE_H
