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


// PassMess includes
#include "libmesh/passmess_assert.h"

// libMesh includes
#include "libmesh/libmesh_common.h"

// C/C++ includes
#ifdef LIBMESH_HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef LIBMESH_HAVE_UNISTD_H
#include <unistd.h>  // needed for getpid()
#endif

#ifdef LIBMESH_HAVE_CSIGNAL
#  include <csignal>
#endif

namespace PassMess
{

void report_here(const char * file, int line, const char * date, const char * time)
{
  std::cerr << "[" << static_cast<std::size_t>(libMesh::global_processor_id()) << "] "
            << file
            << ", line " << line
            << ", compiled " << date
            << " at " << time
            << std::endl;
}


void report_error(const char * file, int line, const char * date, const char * time)
{
  // It is possible to have an error *inside* report_error; e.g. when
  // we start using a PassMess::print_trace.  Don't infinitely recurse.
  static bool reporting_error = false;
  if (reporting_error)
    {
      // I heard you like error reporting, so we put an error report
      // in report_error() so you can report errors from the report.
      std::cerr << "PassMess encountered an error while attempting to report_error." << std::endl;
      return;
    }
  reporting_error = true;

  report_here(file, line, date, time);

  reporting_error = false;
}

} // namespace PassMess
