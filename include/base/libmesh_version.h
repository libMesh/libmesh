// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_LIBMESH_VERSION_H
#define LIBMESH_LIBMESH_VERSION_H

#include "libmesh_config.h"

// C++ includes
#include <string>

// You can use this macro to guard pieces of your application code
// against use in incorrect versions of libmesh.  For example:
// #if LIBMESH_VERSION_LESS_THAN(1,0,0)
// ...
// #elif LIBMESH_VERSION_LESS_THAN(1,1,0)
// ...
// #endif
#define LIBMESH_VERSION_LESS_THAN(major,minor,micro)                    \
  ((LIBMESH_MAJOR_VERSION < (major) ||                                  \
    (LIBMESH_MAJOR_VERSION == (major) && (LIBMESH_MINOR_VERSION < (minor) || \
                                          (LIBMESH_MINOR_VERSION == (minor) && \
                                           LIBMESH_MICRO_VERSION < (micro))))) ? 1 : 0)

namespace libMesh
{
void        libmesh_version_stdout();
int         get_libmesh_version();

/**
 * Specifier for I/O file compatibility features.
 * This only needs to be changed when new restart file
 * functionality is added.
 */
std::string get_io_compatibility_version();
}


#endif // LIBMESH_LIBMESH_VERSION_H
