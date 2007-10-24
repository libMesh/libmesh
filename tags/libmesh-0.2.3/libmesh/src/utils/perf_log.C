// $Id: perf_log.C,v 1.5 2003-01-24 17:24:45 jwpeterson Exp $

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



// C++ includes
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <pwd.h>
#include <sstream>

// Local includes
#include "perf_log.h"




// ------------------------------------------------------------
// PerfLog class member funcions

bool PerfLog::called = false;


PerfLog::PerfLog(std::string cn,
		 const bool le) :
  class_name(cn),
  log_events(le)
{
  if (log_events)
    clear();
};



PerfLog::~PerfLog()
{
  if (log_events)
    print_log();
};



void PerfLog::clear()
{
  if (log_events)
    {
      //  check that all events are closed
      for (std::map<std::string, PerfData>::iterator
	     pos = log.begin(); pos != log.end(); ++pos)
	if (pos->second.open)
	  {
	    std::cout
	      << "ERROR clearning performance log for class "
	      << class_name << std::endl
	      << "event " << pos->first << " is still being monitored!"
	      << std::endl;

	    error();
	  }

      
      gettimeofday (&tstart, NULL);
  
      log.clear();
    }
};


std::string PerfLog::get_info_header() const
{
  std::ostringstream out;
  
  if (log_events)
    {

#ifdef HAVE_LOCALE
      std::ostringstream  dateStr;
      time_t tm         = time(NULL);
      struct tm* tmb    = localtime(&tm);
      std::locale loc;
      TimeIter            begin(dateStr);
      const TimePut& tp = std::use_facet<TimePut>(loc);
      tp.put(begin,
	     dateStr,
	     dateStr.fill(),
	     tmb,
	     'c');
#endif
      
      // Get system information
      struct utsname sysInfo;
      uname(&sysInfo);
      
      // Get user information
      struct passwd* p = getpwuid(getuid());
      out << std::endl
	  << " ---------------------------------------------------------------------" << std::endl
#ifdef HAVE_LOCALE
	  << "| Time:           " << dateStr.str()    << std::endl
#endif
	  << "| OS:             " << sysInfo.sysname  << std::endl
	  << "| HostName:       " << sysInfo.nodename << std::endl
	  << "| OS Release      " << sysInfo.release  << std::endl
	  << "| OS Version:     " << sysInfo.version  << std::endl
	  << "| Machine:        " << sysInfo.machine  << std::endl
	  << "| Username:       " << p->pw_name       << std::endl 
	  << " ---------------------------------------------------------------------" << std::endl;
    }

  return out.str();
}




std::string PerfLog::get_perf_info() const
{
  std::ostringstream out;

#ifndef BROKEN_IOSTREAM
  
  if ((log_events) && (!log.empty()))
    {
      struct timeval tstop;
      
      gettimeofday (&tstop, NULL);
	  
      const double elapsed_time = ((double) (tstop.tv_sec  - tstart.tv_sec)) +
	((double) (tstop.tv_usec - tstart.tv_usec))/1000000.;      
      
      out << " ---------------------------------------------------------------------"  << std::endl;
      out << "| " << class_name << " Performance: elapsed_time=" << elapsed_time       << std::endl;
      out << " ---------------------------------------------------------------------"  << std::endl;
      out << "| ";
      out.width(24);
      out << std::left << "Event";
      
      out.width(8);
      out << std::left << "nCalls";
      
      out.width(12);
      out << std::left << "Total";
      
      out.width(12);
      out << std::left << "Avg";
      
      out.width(12);
      out << std::left << "Percent of";
      
      out << "|" << std::endl;
      
      out << "| ";
      out.width(24);
      out << std::left << "";
      
      out.width(8);
      out << std::left << "";
      
      out.width(12);
      out << std::left << "Time";
      
      out.width(12);
      out << std::left << "Time";
      
      out.width(12);
      out << std::left << "Total Time";
      
      out << "|" << std::endl;
      out << "|---------------------------------------------------------------------|" << std::endl
	  << "|                                                                     |" << std::endl;
      
      
      for (std::map<std::string, PerfData>::const_iterator
	     pos = log.begin(); pos != log.end(); ++pos)
	{
	  const PerfData perf_data = pos->second;

	  // Only print the event if the count is non-zero.
	  if (perf_data.count != 0)
	    {
	      out << "| ";
	      out.width(24);
	      out << std::left << pos->first;
	      
	      out.width(8);
	      out << perf_data.count;
	      
	      out.setf(std::ios::fixed);
	      out.width(12);
	      out.precision(4);
	      out << perf_data.tot_time;
	      
	      out.width(12);
	      out.precision(4);
	      out << perf_data.tot_time/static_cast<double>(perf_data.count);
	      
	      out.width(12);
	      out.precision(2);
	      out << perf_data.tot_time/elapsed_time*100.;
	      
	      out << "|";
	      out << std::endl;
	    }
	}
      
      out << " ---------------------------------------------------------------------" << std::endl;
    }

#endif
  
  return out.str();
}



std::string PerfLog::get_log() const
{
  std::ostringstream out;
  
#ifndef BROKEN_IOSTREAM
  
  if (log_events)
    {
      // Possibly print machine info      
      if (!called)
	{
	  called = true;
	  out << get_info_header();
	}

     
      // Print the log
      if (!log.empty())
	{
	  out << get_perf_info();
	}
    }

#endif
  
  return out.str();
};



void PerfLog::print_log() const
{
  if (log_events)
    std::cout << get_log() << std::endl;
};
