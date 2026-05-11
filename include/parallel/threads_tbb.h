// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_THREADS_TBB_H
#define LIBMESH_THREADS_TBB_H

// Do not try to #include this header directly, it is designed to be
// #included directly by threads.h
#ifndef LIBMESH_SQUASH_HEADER_WARNING
# warning "This file is designed to be included through libmesh/threads.h"
#else

#ifdef LIBMESH_HAVE_TBB_API

// Standard library includes needed for oneTBB compatibility wrappers.
// Placed outside the warning-suppressor block since these are standard headers.
#ifdef LIBMESH_HAVE_ONETBB
#  include <atomic>
#  include <memory>
#  include <mutex>
#  include <thread>
#endif

// libMesh includes
#include "libmesh/ignore_warnings.h"

// Threading building blocks includes — era-specific
#ifdef LIBMESH_HAVE_ONETBB
// oneTBB (>= 2021): tbb_stddef.h, task_scheduler_init.h, atomic.h,
// tbb_thread.h, and recursive_mutex.h no longer exist.
#  include "tbb/version.h"
#  include "tbb/global_control.h"
#else
// Legacy Intel TBB (< 2021)
#  include "tbb/tbb_stddef.h"
#  include "tbb/task_scheduler_init.h"
#  include "tbb/atomic.h"
#  include "tbb/tbb_thread.h"
#  include "tbb/recursive_mutex.h"
#endif

// Headers present in both eras (oneTBB provides tbb/*.h compatibility wrappers)
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/partitioner.h"
#include "tbb/spin_mutex.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_arena.h"

#include "libmesh/restore_warnings.h"

#define TBB_VERSION_LESS_THAN(major,minor)                              \
  ((LIBMESH_DETECTED_TBB_VERSION_MAJOR < (major) ||                     \
    (LIBMESH_DETECTED_TBB_VERSION_MAJOR == (major) && (LIBMESH_DETECTED_TBB_VERSION_MINOR < (minor)))) ? 1 : 0)

// Thread-Local-Storage macros
#define LIBMESH_TLS_TYPE(type)  tbb::enumerable_thread_specific<type>
#define LIBMESH_TLS_REF(value)  (value).local()

namespace libMesh
{

namespace Threads
{

/**
 * Thread object abstraction that provides a basic constructor and
 * support for join().
 *
 * Legacy TBB: tbb::tbb_thread  |  oneTBB: std::thread
 */
#ifndef LIBMESH_HAVE_ONETBB
typedef tbb::tbb_thread Thread;
#else
typedef std::thread Thread;
#endif

/**
 * Scheduler to manage the TBB thread pool.
 *
 * Legacy TBB: thin typedef around tbb::task_scheduler_init.
 * oneTBB: task_scheduler_init was removed; this wrapper uses
 * tbb::global_control to set the same global parallelism limit.
 */
#ifndef LIBMESH_HAVE_ONETBB
typedef tbb::task_scheduler_init task_scheduler_init;
#else
class task_scheduler_init
{
public:
  static const int automatic = -1;

  explicit task_scheduler_init (int n = automatic)
  {
    if (n != automatic && n > 0)
      _gc = std::make_unique<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(n));
  }

  void initialize (int n = automatic)
  {
    if (n != automatic && n > 0)
      _gc = std::make_unique<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(n));
  }

  void terminate () { _gc.reset(); }

private:
  std::unique_ptr<tbb::global_control> _gc;
};
#endif // LIBMESH_HAVE_ONETBB

/**
 * Dummy "splitting object" used to distinguish splitting constructors
 * from copy constructors.
 */
typedef tbb::split split;

/**
 * Execute the provided function object in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_for (const Range & range, const Body & body,
                   unsigned int n_threads = libMesh::n_threads())
{
  libmesh_error_msg_if(n_threads > libMesh::n_threads(),
                       "Requested n_threads (" << n_threads << ") exceeds the "
                       "global thread count (" << libMesh::n_threads() << ").");
  BoolAcquire b(in_threads);

  if (n_threads > 1)
    {
      Threads::RAIIAcquire<int> a(Threads::active_threads, n_threads);

      DisablePerfLogInScope disable_perf;
      if (n_threads == libMesh::n_threads())
        tbb::parallel_for (range, body, tbb::auto_partitioner());
      else
        {
          tbb::task_arena arena(static_cast<int>(n_threads));
          arena.execute([&]{ tbb::parallel_for(range, body, tbb::auto_partitioner()); });
        }
    }
  else
    body(range);
}



/**
 * Execute the provided function object in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_for (const Range & range, const Body & body, const Partitioner & partitioner,
                   unsigned int n_threads = libMesh::n_threads())
{
  libmesh_error_msg_if(n_threads > libMesh::n_threads(),
                       "Requested n_threads (" << n_threads << ") exceeds the "
                       "global thread count (" << libMesh::n_threads() << ").");
  BoolAcquire b(in_threads);

  if (n_threads > 1)
    {
      Threads::RAIIAcquire<int> a(Threads::active_threads, n_threads);

      DisablePerfLogInScope disable_perf;
      if (n_threads == libMesh::n_threads())
        tbb::parallel_for (range, body, partitioner);
      else
        {
          tbb::task_arena arena(static_cast<int>(n_threads));
          arena.execute([&]{ tbb::parallel_for(range, body, partitioner); });
        }
    }
  else
    body(range);
}



/**
 * Execute the provided reduction operation in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_reduce (const Range & range, Body & body,
                      unsigned int n_threads = libMesh::n_threads())
{
  libmesh_error_msg_if(n_threads > libMesh::n_threads(),
                       "Requested n_threads (" << n_threads << ") exceeds the "
                       "global thread count (" << libMesh::n_threads() << ").");
  BoolAcquire b(in_threads);

  if (n_threads > 1)
    {
      Threads::RAIIAcquire<int> a(Threads::active_threads, n_threads);

      DisablePerfLogInScope disable_perf;
      if (n_threads == libMesh::n_threads())
        tbb::parallel_reduce (range, body, tbb::auto_partitioner());
      else
        {
          tbb::task_arena arena(static_cast<int>(n_threads));
          arena.execute([&]{ tbb::parallel_reduce(range, body, tbb::auto_partitioner()); });
        }
    }
  else
    body(range);
}



/**
 * Execute the provided reduction operation in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_reduce (const Range & range, Body & body, const Partitioner & partitioner,
                      unsigned int n_threads = libMesh::n_threads())
{
  libmesh_error_msg_if(n_threads > libMesh::n_threads(),
                       "Requested n_threads (" << n_threads << ") exceeds the "
                       "global thread count (" << libMesh::n_threads() << ").");
  BoolAcquire b(in_threads);

  if (n_threads > 1)
    {
      Threads::RAIIAcquire<int> a(Threads::active_threads, n_threads);

      DisablePerfLogInScope disable_perf;
      if (n_threads == libMesh::n_threads())
        tbb::parallel_reduce (range, body, partitioner);
      else
        {
          tbb::task_arena arena(static_cast<int>(n_threads));
          arena.execute([&]{ tbb::parallel_reduce(range, body, partitioner); });
        }
    }
  else
    body(range);
}



/**
 * Spin mutex.  Implements mutual exclusion by busy-waiting in user
 * space for the lock to be acquired.
 *
 * tbb::spin_mutex (and its scoped_lock) is present in both legacy and oneTBB.
 */
typedef tbb::spin_mutex spin_mutex;

/**
 * Recursive mutex.  Implements mutual exclusion allowing the same thread
 * to acquire the lock multiple times.
 *
 * Legacy TBB: thin typedef around tbb::recursive_mutex.
 * oneTBB: tbb::recursive_mutex was removed; this wrapper provides the same
 * scoped_lock interface using std::recursive_mutex underneath.
 */
#ifndef LIBMESH_HAVE_ONETBB
typedef tbb::recursive_mutex recursive_mutex;
#else
class recursive_mutex
{
public:
  void lock ()   { _m.lock(); }
  void unlock () { _m.unlock(); }

  class scoped_lock
  {
  public:
    scoped_lock () : _rm(nullptr) {}
    explicit scoped_lock (recursive_mutex & rm) : _rm(nullptr) { acquire(rm); }
    ~scoped_lock () { release(); }

    void acquire (recursive_mutex & rm) { _rm = &rm; _rm->lock(); }
    void release () { if (_rm) { _rm->unlock(); _rm = nullptr; } }

  private:
    recursive_mutex * _rm;
  };

private:
  std::recursive_mutex _m;
};
#endif // LIBMESH_HAVE_ONETBB

/**
 * Defines atomic operations which can only be executed on a
 * single thread at a time.  This is used in reference counting,
 * for example, to allow count++/count-- to work.
 *
 * Legacy TBB: inherits tbb::atomic<T>.
 * oneTBB: tbb::atomic was removed; inherits std::atomic<T> instead.
 */
#ifndef LIBMESH_HAVE_ONETBB
template <typename T>
class atomic : public tbb::atomic<T> {};
#else
template <typename T>
class atomic : public std::atomic<T>
{
public:
  atomic () : std::atomic<T>(0) {}
};
#endif // LIBMESH_HAVE_ONETBB

} // namespace Threads

} // namespace libMesh

#endif // LIBMESH_HAVE_TBB_API

#endif // LIBMESH_SQUASH_HEADER_WARNING

#endif // LIBMESH_THREADS_TBB_H
