// $Id: utility.C,v 1.5 2003-01-24 17:24:46 jwpeterson Exp $

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



// System includes
#include <sstream>
#include <iostream>
#include <sys/time.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/utsname.h>


// Local includes
#include "utility.h"



// The system_info function duplicates some of the
// functionality found in the perf_log function.
// This way you can get information about a user's
// system without creating a perf_log object.
std::string Utility::system_info()
{
  std::ostringstream out;
  
#ifdef HAVE_LOCALE
#ifndef BROKEN_IOSTREAM
    
  std::locale loc;
  std::ostringstream  dateStr;
  std::ostreambuf_iterator<char, std::char_traits<char> > begin(dateStr);
  time_t tm                  = time(NULL);
  struct tm* tmb             = localtime(&tm);
  const Utility::TimePut& tp = std::use_facet<TimePut>(loc);
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
  out << std::endl
      << " ---------------------------------------------------------------------" << std::endl
      << "| Time:           " << dateStr.str()    << std::endl
      << "| OS:             " << sysInfo.sysname  << std::endl
      << "| HostName:       " << sysInfo.nodename << std::endl
      << "| OS Release      " << sysInfo.release  << std::endl
      << "| OS Version:     " << sysInfo.version  << std::endl
      << "| Machine:        " << sysInfo.machine  << std::endl
      << "| Username:       " << p->pw_name       << std::endl 
      << " ---------------------------------------------------------------------" << std::endl;

#endif
#endif
  
  return out.str();
}
