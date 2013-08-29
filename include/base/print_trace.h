// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PRINT_TRACE_H
#define LIBMESH_PRINT_TRACE_H

// Local includes
#include "libmesh/libmesh_config.h"

// C++ includes
#include <iostream>

namespace libMesh
{

/*
 * Print a stack trace (for code compiled with gcc)
 */
void print_trace(std::ostream &out = std::cerr);

/*
 * Mostly system independent demangler
 */
std::string demangle(const char *name);

  /**
   * Writes a stack trace to a uniquely named file if
   * --enable-tracefiles has been set by configure, otherwise does
   * nothing. Note that we append to the trace file rather than
   * overwriting it.  This allows multiple traces to be written to the
   * same file.
   */
  void write_traceout();

} // namespace libMesh

#endif // LIBMESH_PRINT_TRACE_H
