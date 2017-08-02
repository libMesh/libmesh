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


// library configuration
#include "libmesh/libmesh_config.h"

// System includes
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>

#ifdef LIBMESH_HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>
#endif

#ifdef LIBMESH_HAVE_PWD_H
#include <pwd.h>
#endif

#ifdef LIBMESH_HAVE_DIRECT_H
#include <direct.h>
#endif

// Local includes
#include "libmesh/utility.h"
#include "libmesh/timestamp.h"

namespace libMesh
{


//-----------------------------------------------------------------------
// Utility members


// The system_info function duplicates some of the
// functionality found in the perf_log function.
// This way you can get information about a user's
// system without creating a perf_log object.
std::string Utility::system_info()
{
  std::ostringstream oss;

  std::string date = Utility::get_timestamp();

#ifdef LIBMESH_HAVE_SYS_UTSNAME_H
  // Get system information
  struct utsname sysInfo;
  uname(&sysInfo);
#endif

  // Get user information
#ifdef LIBMESH_HAVE_GETPWUID
  struct passwd * p = getpwuid(getuid());
#endif


  oss << '\n'
      << " ---------------------------------------------------------------------\n"
      << "| Time:           " << date             << '\n'
#ifdef LIBMESH_HAVE_SYS_UTSNAME_H
      << "| OS:             " << sysInfo.sysname  << '\n'
      << "| HostName:       " << sysInfo.nodename << '\n'
      << "| OS Release      " << sysInfo.release  << '\n'
      << "| OS Version:     " << sysInfo.version  << '\n'
      << "| Machine:        " << sysInfo.machine  << '\n'
#else
      << "| OS:             " << "Unknown"        << '\n'
      << "| HostName:       " << "Unknown"        << '\n'
      << "| OS Release      " << "Unknown"        << '\n'
      << "| OS Version:     " << "Unknown"        << '\n'
      << "| Machine:        " << "Unknown"        << '\n'
#endif
#ifdef LIBMESH_HAVE_GETPWUID
      << "| Username:       " << p->pw_name       << '\n'
#else
      << "| Username:       " << "Unknown"        << '\n'
#endif
      << " ---------------------------------------------------------------------\n";

  return oss.str();
}



#ifdef LIBMESH_USE_COMPLEX_NUMBERS

std::string Utility::complex_filename (const std::string & basename,
                                       const unsigned int r_o_c)
{
  std::string name(basename);

  if (r_o_c == 0)
    name.append(".real");

  else
    name.append(".imag");

  return name;
}



void Utility::prepare_complex_data(const std::vector<Complex> & source,
                                   std::vector<Real> & real_part,
                                   std::vector<Real> & imag_part)
{
  const unsigned int len = source.size();

  real_part.resize(len);
  imag_part.resize(len);

  for (unsigned int i=0; i<len; i++)
    {
      real_part[i] = source[i].real();
      imag_part[i] = source[i].imag();
    }
}

#endif // #ifdef LIBMESH_USE_COMPLEX_NUMBERS


int Utility::mkdir(const char* pathname) {
#if defined(LIBMESH_HAVE_MKDIR)
  return ::mkdir(pathname, 0755);
#elif LIBMESH_HAVE_DECL__MKDIR
  return _mkdir(pathname);
#else
  libmesh_error_msg("Function mkdir not available on this system.");
#endif

}

} // namespace libMesh
