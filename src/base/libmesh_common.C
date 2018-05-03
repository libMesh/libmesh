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


// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/print_trace.h"

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

namespace libMesh
{

namespace MacroFunctions
{
void here(const char * file, int line, const char * date, const char * time)
{
  libMesh::err << "[" << static_cast<std::size_t>(libMesh::global_processor_id()) << "] "
               << file
               << ", line " << line
               << ", compiled " << date
               << " at " << time
               << std::endl;
}



void stop(const char * file, int line, const char * date, const char * time)
{
  if (libMesh::global_n_processors() == 1)
    {
      libMesh::MacroFunctions::here(file, line, date, time);
#if defined(LIBMESH_HAVE_CSIGNAL) && defined(SIGSTOP)
      libMesh::out << "Stopping process " << getpid() << "..." << std::endl;
      std::raise(SIGSTOP);
      libMesh::out << "Continuing process " << getpid() << "..." << std::endl;
#else
      libMesh::out << "WARNING:  libmesh_stop() does not work; no operating system support." << std::endl;
#endif
    }
}


void report_error(const char * file, int line, const char * date, const char * time)
{
  // It is possible to have an error *inside* report_error; e.g. from
  // print_trace.  We don't want to infinitely recurse.
  static bool reporting_error = false;
  if (reporting_error)
    {
      // I heard you like error reporting, so we put an error report
      // in report_error() so you can report errors from the report.
      libMesh::err << "libMesh encountered an error while attempting to report_error." << std::endl;
      return;
    }
  reporting_error = true;

  if (libMesh::global_n_processors() == 1 ||
      libMesh::on_command_line("--print-trace"))
    libMesh::print_trace();
  else
    libMesh::write_traceout();
  libMesh::MacroFunctions::here(file, line, date, time);

  reporting_error = false;
}

} // namespace MacroFunctions
} // namespace libMesh
