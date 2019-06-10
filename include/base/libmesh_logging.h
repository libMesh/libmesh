// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LIBMESH_LOGGING_H
#define LIBMESH_LIBMESH_LOGGING_H


// The library configuration options
#include "libmesh/libmesh_common.h"

#include "libmesh/perf_log.h"

// Two-level macro substitution trick, used to construct a unique
// variable name for a given line.
#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

namespace libMesh
{


// Forward declaration, required when included
// in perf_log.{C,h} because the preceding
// #include "libmesh/perf_log.h" is ineffective.
// Multiple inclusion avoidance problem...
// LIBMESH_PERF_LOG_H already #define'd, but the
// class has not been declared yet!.
class PerfLog;

/**
 * A \p PerfLog object to log performance.  If the library is configured
 * with \p --enable-perflog then it will log key functions.
 */
extern PerfLog perflog;


/**
 * Used for logging something that naturally lasts as long as some
 * enclosing scope, such as the current function.  Makes it very easy
 * to handle multiple return scenarios, since the event is popped in
 * the destructor.  Should not be used directly, instead use the
 * LOG_SCOPE macro, which resolves to nothing at compile time if
 * logging is disabled.
 *
 * \author John Peterson
 * \date 2016
 */
struct PerfItem
{
  PerfItem(const char * label,
           const char * header,
           bool enabled=true,
           PerfLog * my_perflog=&perflog) :
    _label(label),
    _header(header),
    _enabled(enabled),
    _perflog(*my_perflog)
  {
    if (_enabled)
      _perflog.fast_push(label, header);
  }

  ~PerfItem()
  {
    if (_enabled)
      _perflog.fast_pop(_label, _header);
  }

private:
  const char * _label;
  const char * _header;
  bool _enabled;
  PerfLog & _perflog;
};



} // namespace libMesh



// Macros for performance logging.  This allows us
// to add performance monitors to the code without
// impacting performance when performance logging
// is disabled.
#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING

#  define START_LOG(a,b)   { libMesh::perflog.push(a,b); }
#  define STOP_LOG(a,b)    { libMesh::perflog.pop(a,b); }
#ifdef LIBMESH_ENABLE_DEPRECATED
#  define PALIBMESH_USE_LOG(a,b)   { libmesh_deprecated(); }
#  define RESTART_LOG(a,b) { libmesh_deprecated(); }
#endif
#  define LOG_SCOPE(a,b)   libMesh::PerfItem TOKENPASTE2(perf_item_, __LINE__)(a,b);
#  define LOG_SCOPE_IF(a,b,enabled)   libMesh::PerfItem TOKENPASTE2(perf_item_, __LINE__)(a,b,enabled);
#  define LOG_SCOPE_WITH(a,b,logger)   libMesh::PerfItem TOKENPASTE2(perf_item_, __LINE__)(a,b,true,&logger);

#else

#  define START_LOG(a,b)   {}
#  define STOP_LOG(a,b)    {}
#  define PALIBMESH_USE_LOG(a,b)   {}
#  define RESTART_LOG(a,b) {}
#  define LOG_SCOPE(a,b)   {}
#  define LOG_SCOPE_IF(a,b,enabled) {}
#  define LOG_SCOPE_WITH(a,b,logger) {}

#endif





#endif // LIBMESH_LIBMESH_LOGGING_H
