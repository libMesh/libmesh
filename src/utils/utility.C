// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_logging.h"

// System includes
#ifdef LIBMESH_HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef LIBMESH_HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef LIBMESH_HAVE_UNISTD_H
#include <unistd.h> // for getuid(), getpid()
#endif
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


std::string Utility::unzip_file (const std::string & name)
{
  std::ostringstream pid_suffix;
  pid_suffix << '_' << getpid();

  std::string new_name = name;
  if (name.size() - name.rfind(".bz2") == 4)
    {
#ifdef LIBMESH_HAVE_BZIP
      new_name.erase(new_name.end() - 4, new_name.end());
      new_name += pid_suffix.str();
      LOG_SCOPE("system(bunzip2)", "Utility");
      std::string system_string = "bunzip2 -f -k -c ";
      system_string += name + " > " + new_name;
      if (std::system(system_string.c_str()))
        libmesh_file_error(system_string);
#else
      libmesh_error_msg("ERROR: need bzip2/bunzip2 to open .bz2 file " << name);
#endif
    }
  else if (name.size() - name.rfind(".xz") == 3)
    {
#ifdef LIBMESH_HAVE_XZ
      new_name.erase(new_name.end() - 3, new_name.end());
      new_name += pid_suffix.str();
      LOG_SCOPE("system(xz -d)", "Utility");
      std::string system_string = "xz -f -d -k -c ";
      system_string += name + " > " + new_name;
      if (std::system(system_string.c_str()))
        libmesh_file_error(system_string);
#else
      libmesh_error_msg("ERROR: need xz to open .xz file " << name);
#endif
    }
  return new_name;
}




} // namespace libMesh
