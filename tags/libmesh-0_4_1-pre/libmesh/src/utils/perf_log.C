// $Id: perf_log.C,v 1.17 2003-09-02 18:02:45 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "perf_log.h"
#include "o_string_stream.h"




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
  OStringStream out;
  
  if (log_events)
    {

#ifdef HAVE_LOCALE
      OStringStream  dateStr;
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
  OStringStream out;
  
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
      OSSStringleft(out,30,"Event");      
      OSSStringleft(out,8,"nCalls");      
      OSSStringleft(out,12,"Total");      
      OSSStringleft(out,12,"Avg");      
      OSSStringleft(out,13,"Percent of");      
      out << "|" << std::endl;      
      out << "| ";
      OSSStringleft(out,30,"");
      OSSStringleft(out,8,"");
      OSSStringleft(out,12,"Time");
      OSSStringleft(out,12,"Time");
      OSSStringleft(out,13,"Active Time");
      
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
		  OSSStringleft(out,30,pos->first.second);
		}
	      else
		{
		  if (last_header != pos->first.first)
		    {
		      last_header = pos->first.first;

		      out << "| ";
		      OSSStringright(out,76,"|");
		      out << std::endl;
		      
		      out << "| ";
		      OSSStringleft(out,75,pos->first.first);
		      out << "|";
		      out << std::endl;
		    }

		  out << "|   ";
		  OSSStringleft(out,28,pos->first.second);
		}
	      

	      // Print the number of calls to the event
	      OSSInt(out,8,perf_count);

	      // Print the total time spent in the event
	      out.setf(std::ios::fixed);
	      OSSRealleft(out,12,4,perf_time);

	      // Print the average time per function call
	      OSSRealleft(out,12,6,perf_avg_time);

	      // Print the percentage of the time spent in the event
	      OSSRealleft(out,13,2,perf_percent);
	      
	      out << "|";
	      out << std::endl;
	    }
	}
      
      out << " ----------------------------------------------------------------------------" << std::endl;
      out << "| ";
      OSSStringleft(out,30,"Totals:");

      // Print the total number of logged function calls
      OSSInt(out,8,summed_function_calls);

      // Print the total time spent in logged function calls
      out.setf(std::ios::fixed);
      OSSRealleft(out,12,4,summed_total_time);

      // Null, the average time doesn't make sense as a total
      out.width(12);
      out << "";

      // Print the total percentage
      OSSRealleft(out,13,2,summed_percentage);
      
      out << "|";
      out << std::endl;
      
      out << " ----------------------------------------------------------------------------" << std::endl;
    }

  return out.str();
}



std::string PerfLog::get_log() const
{
  OStringStream out;
  
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
  
  return out.str();
}



void PerfLog::print_log() const
{
  if (log_events)
    std::cout << get_log() << std::endl;
}
