// $Id: mesh_logging.h,v 1.2 2003-05-15 23:34:34 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_logging_h__
#define __mesh_logging_h__


// The library configuration options
#include "mesh_common.h"

// Macros for performance logging.  This allows us
// to add performance monitors to the code without
// impacting performance when performance logging
// is disabled.
#ifdef ENABLE_PERFORMANCE_LOGGING

// Note the log is in libMesh, so we need to include it.
#  include "libmesh.h"
#  define START_LOG(a,b)   { libMesh::log.start_event(a,b); }
#  define STOP_LOG(a,b)    { libMesh::log.stop_event(a,b); }
#  define PAUSE_LOG(a,b)   { libMesh::log.pause_event(a,b); }
#  define RESTART_LOG(a,b) { libMesh::log.restart_event(a,b); }

#else

#  define START_LOG(a,b)   {}
#  define STOP_LOG(a,b)    {}
#  define PAUSE_LOG(a,b)   {}
#  define RESTART_LOG(a,b) {}

#endif





#endif // #define __mesh_logging_h__
