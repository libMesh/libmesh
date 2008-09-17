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



// C++ includes
#include <iostream>
#include <iomanip>
#include <ctime>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <pwd.h>
#include <vector>

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
    this->clear();
}



PerfLog::~PerfLog()
{
  if (log_events)
    this->print_log();
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

	    libmesh_error();
	  }

      
      gettimeofday (&tstart, NULL);
  
      log.clear();

      while (!log_stack.empty())
	log_stack.pop();
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
      out << "\n";

      // Construct string stream objects for each of the outputs
      OStringStream pid_stream;
      OStringStream nprocs_stream;
      OStringStream time_stream;
      OStringStream os_stream;
      OStringStream host_stream;
      OStringStream osrel_stream;
      OStringStream osver_stream;
      OStringStream machine_stream;
      OStringStream user_stream;
      OStringStream config_stream;
      

      // Put pointers to these streams in a vector
      std::vector<OStringStream*> v(10);
      v[0] = &pid_stream;
      v[1] = &nprocs_stream;
      v[2] = &time_stream;
      v[3] = &os_stream;
      v[4] = &host_stream;
      v[5] = &osrel_stream;
      v[6] = &osver_stream;
      v[7] = &machine_stream;
      v[8] = &user_stream;
      v[9] = &config_stream;

      // Fill string stream objects
      if (libMesh::n_processors() > 1)
	{
	  pid_stream     << "| Processor id:   " << libMesh::processor_id(); 
	  nprocs_stream  << "| Num Processors: " << libMesh::n_processors();
	}
      
#ifdef HAVE_LOCALE						       
      time_stream    << "| Time:           " << dateStr.str()          ; 
#endif								       
      
      os_stream      << "| OS:             " << sysInfo.sysname        ; 
      host_stream    << "| HostName:       " << sysInfo.nodename       ; 
      osrel_stream   << "| OS Release:     " << sysInfo.release        ; 
      osver_stream   << "| OS Version:     " << sysInfo.version        ; 
      machine_stream << "| Machine:        " << sysInfo.machine        ; 
      user_stream    << "| Username:       " << p->pw_name             ;
      config_stream  << "| Configuration:  " << LIBMESH_CONFIGURE_INFO;
      
      // Find the longest string, use that to set the line length for formatting.
      unsigned int max_length = 0;
      for (unsigned int i=0; i<v.size(); ++i)
	if (v[i]->str().size() > max_length)
	  max_length = v[i]->str().size();
      
      // Print dashed line
      this->_character_line(max_length+2, '-', out);
      out << '\n';

      // Loop over all the strings and print them out with end-formatting
      for (unsigned int i=0; i<v.size(); ++i)
	{
	  if (v[i]->str().size() > 0)
	    {
	      out << v[i]->str();
	      OSSStringright(out, max_length+4 - v[i]->str().size(), "|\n");
	    }
	}

      // Print dashed line
      this->_character_line(max_length+2, '-', out);
      out << '\n';
    }

  return out.str();
}


 

std::string PerfLog::get_perf_info() const
{
  OStringStream out;
  
  if (log_events && !log.empty())
    {
      // Stop timing for this event.
      struct timeval tstop;
      
      gettimeofday (&tstop, NULL);
	  
      const double elapsed_time = (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
				   static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);

      // Figure out the formatting required based on the event names
      // Unsigned ints for each of the column widths
      unsigned int event_col_width            = 30;
      const unsigned int ncalls_col_width     = 10;
      const unsigned int tot_time_col_width   = 12;
      const unsigned int avg_time_col_width   = 12;
      const unsigned int pct_active_col_width = 13;
      
      // Iterator to be used to loop over the map of timed events
      std::map<std::pair<std::string,std::string>, PerfData>::const_iterator pos;
      
      // Reset the event column width based on the longest event name plus
      // a possible 2-character indentation, plus a space.
      for (pos = log.begin(); pos != log.end(); ++pos)
	if (pos->first.second.size()+3 > event_col_width)
	  event_col_width = pos->first.second.size()+3;

      // Set the total width of the column
      const unsigned int total_col_width =
	event_col_width     +
	ncalls_col_width    +
	tot_time_col_width  +
	avg_time_col_width  +
	pct_active_col_width+1;

      // Print dashed line
      out << ' ';
      this->_character_line(total_col_width, '-', out);
      out << '\n';
            
      {
	// Construct temporary message string
	OStringStream temp;
	temp << "| " << label_name << " Performance: Alive time=" << elapsed_time
	     << ", Active time=" << total_time;

	// Get the size of the temporary string
	const unsigned int temp_size = temp.str().size();

	// Send the temporary message to the output
	out << temp.str();

	// If this string is longer than the previously computed total
	// column width, skip the additional formatting... this shouldn't
	// happen often, hopefully.  Add two additional characters for a
	// space and a "|" character at the end.
	if (temp_size < total_col_width+2)
	  {
	    OSSStringright(out, total_col_width-temp_size+2, "|");
	  }
	
	out << '\n';
      }
      
      // Print dashed line
      out << ' ';
      this->_character_line(total_col_width, '-', out);
      out << '\n';

      
      // Write out the header for the events listing
      out << "| ";
      OSSStringleft(out,event_col_width,"Event");      
      OSSStringleft(out,ncalls_col_width,"nCalls");      
      OSSStringleft(out,tot_time_col_width,"Total");      
      OSSStringleft(out,avg_time_col_width,"Avg");      
      OSSStringleft(out,pct_active_col_width,"Percent of");      
      out << "|\n";      
      out << "| ";
      OSSStringleft(out,event_col_width,"");
      OSSStringleft(out,ncalls_col_width,"");
      OSSStringleft(out,tot_time_col_width,"Time");
      OSSStringleft(out,avg_time_col_width,"Time");
      OSSStringleft(out,pct_active_col_width,"Active Time");
      
      out << "|\n|";
      this->_character_line(total_col_width, '-', out);
      out << "|\n|";
      this->_character_line(total_col_width, ' ', out);
      out << "|\n";
      
      unsigned int summed_function_calls = 0;
      double       summed_total_time     = 0;
      double       summed_percentage     = 0;
      
      std::string last_header("");
	  
      for (pos = log.begin(); pos != log.end(); ++pos)
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
		  OSSStringleft(out,event_col_width,pos->first.second);
		}
	      else
		{
		  if (last_header != pos->first.first)
		    {
		      last_header = pos->first.first;

		      // print blank line
		      out << "|";
		      this->_character_line(total_col_width, ' ', out);
		      out << "|\n";

		      // print header name (account for additional space before
		      // the header)
		      out << "| ";
		      OSSStringleft(out, total_col_width-1, pos->first.first);
		      out << "|\n";
		    }

		  out << "|   ";
		  OSSStringleft(out, event_col_width-2, pos->first.second);
		}
	      

	      // Print the number of calls to the event.  
	      OSSInt(out,ncalls_col_width,perf_count);

	      // Print the total time spent in the event
	      out.setf(std::ios::fixed);
	      OSSRealleft(out,tot_time_col_width,4,perf_time);

	      // Print the average time per function call
	      OSSRealleft(out,avg_time_col_width,6,perf_avg_time);

	      // Print the percentage of the time spent in the event
	      OSSRealleft(out,pct_active_col_width,2,perf_percent);
	      
	      out << "|";
	      out << '\n';
	    }
	}

      out << ' ';
      this->_character_line(total_col_width, '-', out);
      out << '\n';
      out << "| ";
      OSSStringleft(out,event_col_width,"Totals:");

      // Print the total number of logged function calls
      // For routines which are called many times, summed_function_calls may
      // exceed 7 digits.  If this happens use, scientific notation.
      if (summed_function_calls < 9999999)
	{
	  OSSInt(out,ncalls_col_width,summed_function_calls);
	}
      
      else
	{
	  out.setf(std::ios::scientific);
	  OSSRealleft(out, ncalls_col_width, 3, static_cast<Real>(summed_function_calls));
	  out.unsetf(std::ios::scientific);
	}
	
      // Print the total time spent in logged function calls
      out.setf(std::ios::fixed);
      OSSRealleft(out,tot_time_col_width,4,summed_total_time);

      // Null, the average time doesn't make sense as a total
      out.width(avg_time_col_width);
      out << "";

      // Print the total percentage
      OSSRealleft(out,pct_active_col_width,2,summed_percentage);
      
      out << "|\n ";
      this->_character_line(total_col_width, '-', out);
      out << '\n';
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
    {
      // Check to see if the log_string is empty, and if so,
      // avoid printing an unnecessary newline.
      std::string log_string = this->get_log();
      if (log_string.size() > 0)
	std::cout << log_string << std::endl;
    }
}



void PerfLog::start_event(const std::string &label,
			  const std::string &header)
{
  this->push(label,header);
}



void PerfLog::stop_event(const std::string &label,
			 const std::string &header)
{
  this->pop(label,header);
}



void PerfLog::pause_event(const std::string &,
			  const std::string &)
{  
  // nothing to do.  pushing the next object on the stack will handle it
}



void PerfLog::restart_event(const std::string &,
			    const std::string &)
{
  // nothing to do.  popping the top off the stack will handle it.
}




void PerfLog::_character_line(const unsigned int n,
			      const char c,
			      OStringStream& out) const
{
  for (unsigned int i=0; i<n; ++i)
    out << c;
}
