// $Id: perf_log.C,v 1.12 2003-03-03 18:03:39 benkirk Exp $

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


PerfLog::PerfLog(const std::string& ln,
		 const bool le) :
  label_name(ln),
  log_events(le),
  total_time(0.)
{
  if (log_events)
    clear();
}



PerfLog::~PerfLog()
{
  if (log_events)
    print_log();
}



void PerfLog::clear()
{
  if (log_events)
    {
      //  check that all events are closed
      for (std::map<std::pair<std::string,std::string>, PerfData>::iterator
	     pos = log.begin(); pos != log.end(); ++pos)
	if (pos->second.open)
	  {
	    std::cout
	      << "ERROR clearning performance log for class "
	      << label_name << std::endl
	      << "event " << pos->first.second << " is still being monitored!"
	      << std::endl;

	    error();
	  }

      
      gettimeofday (&tstart, NULL);
  
      log.clear();
    }
}


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
	  << " ----------------------------------------------------------------------------" << std::endl;
      
      if (libMeshBase::n_processors() > 1)
	{
      out << "| Processor id:   " << libMeshBase::processor_id() << std::endl
	  << "| Num Processors: " << libMeshBase::n_processors() << std::endl;
	}
      
#ifdef HAVE_LOCALE
      out << "| Time:           " << dateStr.str()               << std::endl;
#endif
      out << "| OS:             " << sysInfo.sysname             << std::endl
	  << "| HostName:       " << sysInfo.nodename            << std::endl
	  << "| OS Release      " << sysInfo.release             << std::endl
	  << "| OS Version:     " << sysInfo.version             << std::endl
	  << "| Machine:        " << sysInfo.machine             << std::endl
	  << "| Username:       " << p->pw_name                  << std::endl 
	  << " ----------------------------------------------------------------------------" << std::endl;
    }

  return out.str();
}




std::string PerfLog::get_perf_info() const
{
  std::ostringstream out;

#ifndef BROKEN_IOSTREAM
  
  if (log_events && !log.empty())
    {
      struct timeval tstop;
      
      gettimeofday (&tstop, NULL);
	  
      const double elapsed_time = (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
				   static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);

      out << " ----------------------------------------------------------------------------"  << std::endl;
      out << "| " << label_name << " Performance: Alive time=" << elapsed_time
	  << ", Active time=" << total_time << std::endl;
      out << " ----------------------------------------------------------------------------"  << std::endl;
      out << "| ";
      out.width(30);
      out << std::left << "Event";
      
      out.width(8);
      out << std::left << "nCalls";
      
      out.width(12);
      out << std::left << "Total";
      
      out.width(12);
      out << std::left << "Avg";
      
      out.width(13);
      out << std::left << "Percent of";
      
      out << "|" << std::endl;
      
      out << "| ";
      out.width(30);
      out << std::left << "";
      
      out.width(8);
      out << std::left << "";
      
      out.width(12);
      out << std::left << "Time";
      
      out.width(12);
      out << std::left << "Time";
      
      out.width(13);
      out << std::left << "Active Time";
      
      out << "|" << std::endl;
      out << "|----------------------------------------------------------------------------|" << std::endl
	  << "|                                                                            |" << std::endl;
      
      
      unsigned int summed_function_calls = 0;
      double       summed_total_time     = 0;
      double       summed_percentage     = 0;
      
      std::string last_header("");
	  
      for (std::map<std::pair<std::string,std::string>, PerfData>::const_iterator
	     pos = log.begin(); pos != log.end(); ++pos)
	{
	  const PerfData& perf_data = pos->second;

	  // Only print the event if the count is non-zero.
	  if (perf_data.count != 0)
	    {
	      const unsigned int perf_count    = perf_data.count;
	      const double       perf_time     = perf_data.tot_time;
	      const double       perf_avg_time = perf_time / static_cast<double>(perf_count);
	      const double       perf_percent  = (total_time != 0.) ? perf_time / total_time * 100. : 0.;

	      summed_function_calls += perf_count;
	      summed_total_time     += perf_time;
	      summed_percentage     += perf_percent;

	      // Print the event name
	      if (pos->first.first == "")
		{
		  out << "| ";
		  out.width(30);
		  out << std::left << pos->first.second;
		}
	      else
		{
		  if (last_header != pos->first.first)
		    {
		      last_header = pos->first.first;

		      out << "| ";
		      out.width(76);
		      out << std::right << "|" << std::endl;
		      
		      out << "| ";
		      out.width(75);
		      out << std::left << pos->first.first;
		      out << std::right << "|" << std::endl;
		    }

		  out << "|   ";
		  out.width(28);
		  out << std::left << pos->first.second;
		}
	      

	      // Print the number of calls to the event
	      out.width(8);
	      out << perf_count;

	      // Print the total time spent in the event
	      out.setf(std::ios::fixed);
	      out.width(12);
	      out.precision(4);
	      out << perf_time;

	      // Print the average time per function call
	      out.width(12);
	      out.precision(6);
	      out << perf_avg_time;

	      // Print the percentage of the time spent in the event
	      out.width(13);
	      out.precision(2);
	      out << perf_percent;
	      
	      out << "|";
	      out << std::endl;
	    }
	}
      
      out << " ----------------------------------------------------------------------------" << std::endl;
      out << "| ";
      out.width(30);
      out << std::left << "Totals:";

      // Print the total number of logged function calls
      out.width(8);
      out << summed_function_calls;

      // Print the total time spent in logged function calls
      out.setf(std::ios::fixed);
      out.width(12);
      out.precision(4);
      out << summed_total_time;

      // Null, the average time doesn't make sense as a total
      out.width(12);
      out << "";

      // Print the total percentage
      out.width(13);
      out.precision(2);
      out << summed_percentage;
      
      out << "|";
      out << std::endl;
      
      out << " ----------------------------------------------------------------------------" << std::endl;
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
      // Only print the log
      // if it isn't empty
      if (!log.empty())
	{
	  // Possibly print machine info,
	  // but only do this once
	  if (!called)
	    {
	      called = true;
	      out << get_info_header();
	    }	  
	  out << get_perf_info();
	}
    }

#endif
  
  return out.str();
}



void PerfLog::print_log() const
{
  if (log_events)
    std::cout << get_log() << std::endl;
}
