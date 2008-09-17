// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "libmesh_config.h"

// System includes
#include <sys/time.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/utsname.h>

#ifdef LIBMESH_HAVE_LOCALE
# include <locale>
#endif

// Local includes
#include "o_string_stream.h"
#include "utility.h"



//-----------------------------------------------------------------------
// Utility members


// The system_info function duplicates some of the
// functionality found in the perf_log function.
// This way you can get information about a user's
// system without creating a perf_log object.
std::string Utility::system_info()
{
  OStringStream out;
  
#ifdef LIBMESH_HAVE_LOCALE
    
  typedef std::ostreambuf_iterator<char, std::char_traits<char> > TimeIter;
  typedef std::time_put<char, TimeIter>                           TimePut;
 
  std::locale loc;
  OStringStream dateStr;
  std::ostreambuf_iterator<char, std::char_traits<char> > begin(dateStr);
  time_t tm         = time(NULL);
  struct tm* tmb    = localtime(&tm);
  const TimePut& tp = std::use_facet<TimePut>(loc);
  tp.put(begin,
	 dateStr,
	 dateStr.fill(),
	 tmb,
	 'c');
  
  // Get system information
  struct utsname sysInfo;
  uname(&sysInfo);
  
  // Get user information
  struct passwd* p = getpwuid(getuid());
  out << '\n'
      << " ---------------------------------------------------------------------\n"
      << "| Time:           " << dateStr.str()    << '\n'
      << "| OS:             " << sysInfo.sysname  << '\n'
      << "| HostName:       " << sysInfo.nodename << '\n'
      << "| OS Release      " << sysInfo.release  << '\n'
      << "| OS Version:     " << sysInfo.version  << '\n'
      << "| Machine:        " << sysInfo.machine  << '\n'
      << "| Username:       " << p->pw_name       << '\n' 
      << " ---------------------------------------------------------------------\n";

#endif
  
  return out.str();
}



#ifdef LIBMESH_USE_COMPLEX_NUMBERS

std::string Utility::complex_filename (const std::string& basename,
				       const unsigned int r_o_c)
{
  std::string name(basename);
  
  if (r_o_c == 0)
    name.append(".real");
  
  else
    name.append(".imag");
  
  return name;
}



void Utility::prepare_complex_data(const std::vector<Complex>& source,
				   std::vector<Real>& real_part,
				   std::vector<Real>& imag_part)
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
