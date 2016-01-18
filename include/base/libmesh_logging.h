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



#ifndef LIBMESH_LIBMESH_LOGGING_H
#define LIBMESH_LIBMESH_LOGGING_H


// The library configuration options
#include "libmesh/libmesh_common.h"

#include "libmesh/perf_log.h"

namespace libMesh
{


// Forward declaration, required when included
// in perf_log.{C,h} because the preceeding
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


} // namespace libMesh



// Macros for performance logging.  This allows us
// to add performance monitors to the code without
// impacting performance when performance logging
// is disabled.
#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING

#  define START_LOG(a,b)   { libMesh::perflog.push(a,b); }
#  define STOP_LOG(a,b)    { libMesh::perflog.pop(a,b); }
#  define PALIBMESH_USE_LOG(a,b)   { libmesh_deprecated(); }
#  define RESTART_LOG(a,b) { libmesh_deprecated(); }

#else

#  define START_LOG(a,b)   {}
#  define STOP_LOG(a,b)    {}
#  define PALIBMESH_USE_LOG(a,b)   {}
#  define RESTART_LOG(a,b) {}

#endif





#endif // LIBMESH_LIBMESH_LOGGING_H
