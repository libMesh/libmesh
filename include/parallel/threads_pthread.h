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

#ifndef LIBMESH_THREADS_PTHREAD_H
#define LIBMESH_THREADS_PTHREAD_H

// Do not try to #include this header directly, it is designed to be
// #included directly by threads.h
#ifdef LIBMESH_HAVE_PTHREAD

// C++ includes
#ifdef LIBMESH_HAVE_CXX11_THREAD
# include <thread>
#endif

#include "libmesh/libmesh_logging.h"
#include <pthread.h>
#include <algorithm>
#include <vector>

#ifdef __APPLE__
#include <libkern/OSAtomic.h>
#endif

// Thread-Local-Storage macros
#ifdef LIBMESH_HAVE_CXX11_THREAD
#  define LIBMESH_TLS_TYPE(type)  thread_local type
#  define LIBMESH_TLS_REF(value)  (value)
#else // Maybe support gcc __thread eventually?
#  define LIBMESH_TLS_TYPE(type)  type
#  define LIBMESH_TLS_REF(value)  (value)
#endif

namespace libMesh
{

namespace Threads
{


#ifdef LIBMESH_HAVE_CXX11_THREAD
/**
 * Use std::thread when available.
 */
typedef std::thread Thread;

#else

/**
 * Use the non-concurrent placeholder.
 */
typedef NonConcurrentThread Thread;

#endif // LIBMESH_HAVE_CXX11_THREAD


/**
 * Spin mutex.  Implements mutual exclusion by busy-waiting in user
 * space for the lock to be acquired.
 */
#ifdef __APPLE__
class spin_mutex
{
public:
  spin_mutex() : slock(0) {} // The convention is that the lock being zero is _unlocked_
  ~spin_mutex() {}

  void lock () { OSSpinLockLock(&slock); }
  void unlock () { OSSpinLockUnlock(&slock); }

  class scoped_lock
  {
  public:
    scoped_lock () : smutex(libmesh_nullptr) {}
    explicit scoped_lock ( spin_mutex & in_smutex ) : smutex(&in_smutex) { smutex->lock(); }

    ~scoped_lock () { release(); }

    void acquire ( spin_mutex & in_smutex ) { smutex = &in_smutex; smutex->lock(); }
    void release () { if (smutex) smutex->unlock(); smutex = libmesh_nullptr; }

  private:
    spin_mutex * smutex;
  };

private:
  OSSpinLock slock;
};
#else
class spin_mutex
{
public:
  // Might want to use PTHREAD_MUTEX_ADAPTIVE_NP on Linux, but it's not available on OSX.
  spin_mutex() { pthread_spin_init(&slock, PTHREAD_PROCESS_PRIVATE); }
  ~spin_mutex() { pthread_spin_destroy(&slock); }

  void lock () { pthread_spin_lock(&slock); }
  void unlock () { pthread_spin_unlock(&slock); }

  class scoped_lock
  {
  public:
    scoped_lock () : smutex(libmesh_nullptr) {}
    explicit scoped_lock ( spin_mutex & in_smutex ) : smutex(&in_smutex) { smutex->lock(); }

    ~scoped_lock () { release(); }

    void acquire ( spin_mutex & in_smutex ) { smutex = &in_smutex; smutex->lock(); }
    void release () { if (smutex) smutex->unlock(); smutex = libmesh_nullptr; }

  private:
    spin_mutex * smutex;
  };

private:
  pthread_spinlock_t slock;
};
#endif // __APPLE__



/**
 * Recursive mutex.  Implements mutual exclusion by busy-waiting in user
 * space for the lock to be acquired.
 */
class recursive_mutex
{
public:
  // Might want to use PTHREAD_MUTEX_ADAPTIVE_NP on Linux, but it's not available on OSX.
  recursive_mutex()
  {
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&mutex, &attr);
  }
  ~recursive_mutex() { pthread_mutex_destroy(&mutex); }

  void lock () { pthread_mutex_lock(&mutex); }
  void unlock () { pthread_mutex_unlock(&mutex); }

  class scoped_lock
  {
  public:
    scoped_lock () : rmutex(libmesh_nullptr) {}
    explicit scoped_lock ( recursive_mutex & in_rmutex ) : rmutex(&in_rmutex) { rmutex->lock(); }

    ~scoped_lock () { release(); }

    void acquire ( recursive_mutex & in_rmutex ) { rmutex = &in_rmutex; rmutex->lock(); }
    void release () { if (rmutex) rmutex->unlock(); rmutex = libmesh_nullptr; }

  private:
    recursive_mutex * rmutex;
  };

private:
  pthread_mutex_t mutex;
  pthread_mutexattr_t attr;
};

template <typename Range>
unsigned int num_pthreads(Range & range)
{
  unsigned int min = std::min((std::size_t)libMesh::n_threads(), range.size());
  return min > 0 ? min : 1;
}

template <typename Range, typename Body>
class RangeBody
{
public:
  Range * range;
  Body * body;
};

template <typename Range, typename Body>
void * run_body(void * args)
{
  RangeBody<Range, Body> * range_body = (RangeBody<Range, Body> *)args;

  Body & body = *range_body->body;
  Range & range = *range_body->range;

  body(range);

  return libmesh_nullptr;
}

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
 * Execute the provided function object in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_for (const Range & range, const Body & body)
{
  Threads::BoolAcquire b(Threads::in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
  const bool logging_was_enabled = libMesh::perflog.logging_enabled();

  if (libMesh::n_threads() > 1)
    libMesh::perflog.disable_logging();
#endif

  unsigned int n_threads = num_pthreads(range);

  std::vector<Range *> ranges(n_threads);
  std::vector<RangeBody<const Range, const Body> > range_bodies(n_threads);
  std::vector<pthread_t> threads(n_threads);

  // Create the ranges for each thread
  unsigned int range_size = range.size() / n_threads;

  typename Range::const_iterator current_beginning = range.begin();

  for (unsigned int i=0; i<n_threads; i++)
    {
      unsigned int this_range_size = range_size;

      if (i+1 == n_threads)
        this_range_size += range.size() % n_threads; // Give the last one the remaining work to do

      ranges[i] = new Range(range, current_beginning, current_beginning + this_range_size);

      current_beginning = current_beginning + this_range_size;
    }

  // Create the RangeBody arguments
  for (unsigned int i=0; i<n_threads; i++)
    {
      range_bodies[i].range = ranges[i];
      range_bodies[i].body = &body;
    }

  // Create the threads.  It may seem redundant to wrap a pragma in
  // #ifdefs... but GCC warns about an "unknown pragma" if it
  // encounters this line of code when -fopenmp is not passed to the
  // compiler.
#ifdef LIBMESH_HAVE_OPENMP
#pragma omp parallel for schedule (static)
#endif
  for (unsigned int i=0; i<n_threads; i++)
    {
#if !LIBMESH_HAVE_OPENMP
      pthread_create(&threads[i], libmesh_nullptr, &run_body<Range, Body>, (void *)&range_bodies[i]);
#else
      run_body<Range, Body>((void *)&range_bodies[i]);
#endif
    }

#if !LIBMESH_HAVE_OPENMP
  // Wait for them to finish

  // The use of 'int' instead of unsigned for the iteration variable
  // is deliberate here.  This is an OpenMP loop, and some older
  // compilers warn when you don't use int for the loop index.  The
  // reason has to do with signed vs. unsigned integer overflow
  // behavior and optimization.
  // http://blog.llvm.org/2011/05/what-every-c-programmer-should-know.html
  for (int i=0; i<static_cast<int>(n_threads); i++)
      pthread_join(threads[i], libmesh_nullptr);
#endif

  // Clean up
  for (unsigned int i=0; i<n_threads; i++)
    delete ranges[i];

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
  if (libMesh::n_threads() > 1 && logging_was_enabled)
    libMesh::perflog.enable_logging();
#endif
}

/**
 * Execute the provided function object in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_for (const Range & range, const Body & body, const Partitioner &)
{
  parallel_for (range, body);
}

/**
 * Execute the provided reduction operation in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_reduce (const Range & range, Body & body)
{
  Threads::BoolAcquire b(Threads::in_threads);

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
  const bool logging_was_enabled = libMesh::perflog.logging_enabled();

  if (libMesh::n_threads() > 1)
    libMesh::perflog.disable_logging();
#endif

  unsigned int n_threads = num_pthreads(range);

  std::vector<Range *> ranges(n_threads);
  std::vector<Body *> bodies(n_threads);
  std::vector<RangeBody<Range, Body> > range_bodies(n_threads);

  // Create copies of the body for each thread
  bodies[0] = &body; // Use the original body for the first one
  for (unsigned int i=1; i<n_threads; i++)
    bodies[i] = new Body(body, Threads::split());

  // Create the ranges for each thread
  unsigned int range_size = range.size() / n_threads;

  typename Range::const_iterator current_beginning = range.begin();

  for (unsigned int i=0; i<n_threads; i++)
    {
      unsigned int this_range_size = range_size;

      if (i+1 == n_threads)
        this_range_size += range.size() % n_threads; // Give the last one the remaining work to do

      ranges[i] = new Range(range, current_beginning, current_beginning + this_range_size);

      current_beginning = current_beginning + this_range_size;
    }

  // Create the RangeBody arguments
  for (unsigned int i=0; i<n_threads; i++)
    {
      range_bodies[i].range = ranges[i];
      range_bodies[i].body = bodies[i];
    }

  // Create the threads
  std::vector<pthread_t> threads(n_threads);

  // It may seem redundant to wrap a pragma in #ifdefs... but GCC
  // warns about an "unknown pragma" if it encounters this line of
  // code when -fopenmp is not passed to the compiler.
#ifdef LIBMESH_HAVE_OPENMP
#pragma omp parallel for schedule (static)
#endif
  // The use of 'int' instead of unsigned for the iteration variable
  // is deliberate here.  This is an OpenMP loop, and some older
  // compilers warn when you don't use int for the loop index.  The
  // reason has to do with signed vs. unsigned integer overflow
  // behavior and optimization.
  // http://blog.llvm.org/2011/05/what-every-c-programmer-should-know.html
  for (int i=0; i<static_cast<int>(n_threads); i++)
    {
#if !LIBMESH_HAVE_OPENMP
      pthread_create(&threads[i], libmesh_nullptr, &run_body<Range, Body>, (void *)&range_bodies[i]);
#else
      run_body<Range, Body>((void *)&range_bodies[i]);
#endif
    }

#if !LIBMESH_HAVE_OPENMP
  // Wait for them to finish
  for (unsigned int i=0; i<n_threads; i++)
      pthread_join(threads[i], libmesh_nullptr);
#endif

  // Join them all down to the original Body
  for (unsigned int i=n_threads-1; i != 0; i--)
    bodies[i-1]->join(*bodies[i]);

  // Clean up
  for (unsigned int i=1; i<n_threads; i++)
    delete bodies[i];
  for (unsigned int i=0; i<n_threads; i++)
    delete ranges[i];

#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
  if (libMesh::n_threads() > 1 && logging_was_enabled)
    libMesh::perflog.enable_logging();
#endif
}

/**
 * Execute the provided reduction operation in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_reduce (const Range & range, Body & body, const Partitioner &)
{
  parallel_reduce(range, body);
}


/**
 * Defines atomic operations which can only be executed on a
 * single thread at a time.
 */
template <typename T>
class atomic
{
public:
  atomic () : val(0) {}
  operator T () { return val; }

  T operator=( T value )
  {
    spin_mutex::scoped_lock lock(smutex);
    val = value;
    return val;
  }

  atomic<T> & operator=( const atomic<T> & value )
  {
    spin_mutex::scoped_lock lock(smutex);
    val = value;
    return *this;
  }


  T operator+=(T value)
  {
    spin_mutex::scoped_lock lock(smutex);
    val += value;
    return val;
  }

  T operator-=(T value)
  {
    spin_mutex::scoped_lock lock(smutex);
    val -= value;
    return val;
  }

  T operator++()
  {
    spin_mutex::scoped_lock lock(smutex);
    val++;
    return val;
  }

  T operator++(int)
  {
    spin_mutex::scoped_lock lock(smutex);
    val++;
    return val;
  }

  T operator--()
  {
    spin_mutex::scoped_lock lock(smutex);
    val--;
    return val;
  }

  T operator--(int)
  {
    spin_mutex::scoped_lock lock(smutex);
    val--;
    return val;
  }

private:
  T val;
  spin_mutex smutex;
};

} // namespace Threads

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PTHREAD

#endif // LIBMESH_THREADS_PTHREAD_H
