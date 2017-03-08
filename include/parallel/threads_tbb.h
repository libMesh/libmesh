// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh includes
#include "libmesh/libmesh_logging.h"

#include "libmesh/ignore_warnings.h"

// Threading building blocks includes
#include "tbb/tbb_stddef.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/partitioner.h"
#include "tbb/spin_mutex.h"
#include "tbb/recursive_mutex.h"
#include "tbb/atomic.h"
#include "tbb/tbb_thread.h"
#include "tbb/enumerable_thread_specific.h"

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
 */
typedef tbb::tbb_thread Thread;

/**
 * Scheduler to manage threads.
 */
typedef tbb::task_scheduler_init task_scheduler_init;

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
void parallel_for (const Range & range, const Body & body)
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



/**
 * Execute the provided function object in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_for (const Range & range, const Body & body, const Partitioner & partitioner)
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



/**
 * Execute the provided reduction operation in parallel on the specified
 * range.
 */
template <typename Range, typename Body>
inline
void parallel_reduce (const Range & range, Body & body)
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



/**
 * Execute the provided reduction operation in parallel on the specified
 * range with the specified partitioner.
 */
template <typename Range, typename Body, typename Partitioner>
inline
void parallel_reduce (const Range & range, Body & body, const Partitioner & partitioner)
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



/**
 * Spin mutex.  Implements mutual exclusion by busy-waiting in user
 * space for the lock to be acquired.
 */
typedef tbb::spin_mutex spin_mutex;

/**
 * Recursive mutex.  Implements mutual exclusion by busy-waiting in user
 * space for the lock to be acquired.  The same thread can aquire the
 * same lock multiple times
 */
typedef tbb::recursive_mutex recursive_mutex;

/**
 * Defines atomic operations which can only be executed on a
 * single thread at a time.  This is used in reference counting,
 * for example, to allow count++/count-- to work.
 */
template <typename T>
class atomic : public tbb::atomic<T> {};

} // namespace Threads

} // namespace libMesh

#endif // LIBMESH_HAVE_TBB_API

#endif // LIBMESH_SQUASH_HEADER_WARNING

#endif // LIBMESH_THREADS_TBB_H
