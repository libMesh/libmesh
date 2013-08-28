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


// libmesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/print_trace.h"

// C/C++ includes
#include <unistd.h>  // needed for getpid()

#ifdef LIBMESH_HAVE_CSIGNAL
#  include <csignal>
#endif

namespace libMesh
{

namespace MacroFunctions
{
  void here(const char* file, int line, const char* date, const char* time)
  {
    libMesh::err << "[" << static_cast<std::size_t>(libMesh::processor_id()) << "] "
                 << file
                 << ", line " << line
                 << ", compiled " << date
                 << " at " << time
                 << std::endl;
  }



  void stop(const char* file, int line, const char* date, const char* time)
  {
    if (libMesh::n_processors() == 1)
      {
        libMesh::MacroFunctions::here(file, line, date, time);
#ifdef LIBMESH_HAVE_CSIGNAL
        libMesh::out << "Stopping process " << getpid() << "..." << std::endl;
        std::raise(SIGSTOP);
        libMesh::out << "Continuing process " << getpid() << "..." << std::endl;
#else
        libMesh::out << "WARNING:  libmesh_stop() does not work without the <csignal> header file!" << std::endl;
#endif
      }
  }


  void report_error(const char* file, int line, const char* date, const char* time)
  {
    if (libMesh::n_processors() == 1)
      libMesh::print_trace();
    else
      libMesh::write_traceout();
    libMesh::MacroFunctions::here(file, line, date, time);
  }

} // namespace MacroFunctions
} // namespace libMesh
