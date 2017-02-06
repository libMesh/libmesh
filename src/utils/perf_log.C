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



// C++ includes
#include <iostream>
#include <iomanip>
#include <ctime>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <pwd.h>
#include <vector>
#include <sstream>

// Local includes
#include "libmesh/perf_log.h"
#include "libmesh/timestamp.h"

namespace libMesh
{


// ------------------------------------------------------------
// PerfLog class member funcions

bool PerfLog::called = false;


PerfLog::PerfLog(const std::string & ln,
                 const bool le) :
  label_name(ln),
  log_events(le),
  total_time(0.)
{
  gettimeofday (&tstart, libmesh_nullptr);

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
          libmesh_error_msg("ERROR clearning performance log for class " \
                            << label_name                             \
                            << "\nevent "                             \
                            << pos->first.second                      \
                            << " is still being monitored!");

      gettimeofday (&tstart, libmesh_nullptr);

      log.clear();

      while (!log_stack.empty())
        log_stack.pop();
    }
}


std::string PerfLog::get_info_header() const
{
  std::ostringstream oss;

  if (log_events)
    {
      std::string date = Utility::get_timestamp();

      // Get system information
      struct utsname sysInfo;
      uname(&sysInfo);

      // Get user information
      //
      // Some systems, for example Crays, actually have getpwuid on the head-node
      // but (if I understand correctly) a dynamically-linked glibc is not available
      // on the backend, which is running a reduced operating system like Compute
      // Node Linux.  Thus functions like getpwuid cannot be called.  This makes
      // automatically testing for the existence of getpwuid on the login node
      // difficult.  The configure test would work on the login node but will fail
      // on the backend.  Hence we have added a configure flag, --disable-getpwuid,
      // to manually turn this off.
#ifdef LIBMESH_HAVE_GETPWUID
      struct passwd * p = getpwuid(getuid());
#endif
      oss << "\n";

      // Construct string stream objects for each of the outputs
      std::ostringstream
        pid_stream,
        nprocs_stream,
        time_stream,
        os_stream,
        host_stream,
        osrel_stream,
        osver_stream,
        machine_stream,
        user_stream;


      // Put pointers to these streams in a vector
      std::vector<std::ostringstream*> v;
      v.push_back(&pid_stream);
      v.push_back(&nprocs_stream);
      v.push_back(&time_stream);
      v.push_back(&os_stream);
      v.push_back(&host_stream);
      v.push_back(&osrel_stream);
      v.push_back(&osver_stream);
      v.push_back(&machine_stream);
      v.push_back(&user_stream);

      // Fill string stream objects
      if (libMesh::global_n_processors() > 1)
        {
          pid_stream     << "| Processor id:   " << libMesh::global_processor_id();
          nprocs_stream  << "| Num Processors: " << libMesh::global_n_processors();
        }

      time_stream    << "| Time:           " << date                   ;
      os_stream      << "| OS:             " << sysInfo.sysname        ;
      host_stream    << "| HostName:       " << sysInfo.nodename       ;
      osrel_stream   << "| OS Release:     " << sysInfo.release        ;
      osver_stream   << "| OS Version:     " << sysInfo.version        ;
      machine_stream << "| Machine:        " << sysInfo.machine        ;
      user_stream    << "| Username:       ";
#ifdef LIBMESH_HAVE_GETPWUID
      if (p && p->pw_name)
        user_stream  << p->pw_name;
      else
#endif
        user_stream  << "Unknown";

      // Parse the LIBMESH_CONFIGURE_INFO string literal before using it in PerfLog output
      std::string libmesh_configure_info(LIBMESH_CONFIGURE_INFO);
      std::vector<std::string> parsed_libmesh_configure_info;
      this->split_on_whitespace(libmesh_configure_info,
                                parsed_libmesh_configure_info);

      // There should always be at at least one entry in
      // parsed_libmesh_configure_info, even if the user just ran
      // ../configure.
      libmesh_assert_greater (parsed_libmesh_configure_info.size(), 0);

      // Find the longest string in all the streams
      unsigned int max_length = 0;
      for (std::size_t i=0; i<v.size(); ++i)
        if (v[i]->str().size() > max_length)
          max_length = cast_int<unsigned int>
            (v[i]->str().size());

      // Find the longest string in the parsed_libmesh_configure_info
      for (std::size_t i=0; i<parsed_libmesh_configure_info.size(); ++i)
        if (parsed_libmesh_configure_info[i].size() > max_length)
          max_length = cast_int<unsigned int>
            (parsed_libmesh_configure_info[i].size());

      // Print dashed line for the header
      oss << ' '
          << std::string(max_length+1, '-')
          << '\n';

      // Loop over all the strings and add end formatting
      for (std::size_t i=0; i<v.size(); ++i)
        {
          if (v[i]->str().size())
            oss << v[i]->str()
                << std::setw (cast_int<int>
                              (max_length + 4 - v[i]->str().size()))
                << std::right
                << "|\n";
        }

      // Print out configuration header plus first parsed string.  The
      // magic number 18 below accounts for the length of the word
      // 'Configuration'.
      oss << "| Configuration:  "
          << parsed_libmesh_configure_info[0]
          << std::setw (cast_int<int>
                        (max_length + 4 -
                         parsed_libmesh_configure_info[0].size() - 18))
          << std::right
          << "|\n";

      // Loop over the parsed_libmesh_configure_info and add end formatting.  The magic
      // number 3 below accounts for the leading 'pipe' character and indentation
      for (std::size_t i=1; i<parsed_libmesh_configure_info.size(); ++i)
        {
          oss << "|  "
              << parsed_libmesh_configure_info[i]
              << std::setw (cast_int<int>
                            (max_length + 4 -
                             parsed_libmesh_configure_info[i].size() - 3))
              << std::right
              << "|\n";
        }


      // Print dashed line
      oss << ' '
          << std::string(max_length+1, '-')
          << '\n';
    }

  return oss.str();
}




std::string PerfLog::get_perf_info() const
{
  std::ostringstream oss;

  if (log_events && !log.empty())
    {
      // Stop timing for this event.
      struct timeval tstop;

      gettimeofday (&tstop, libmesh_nullptr);

      const double elapsed_time = (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
                                   static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);

      // Figure out the formatting required based on the event names
      // Unsigned ints for each of the column widths
      unsigned int event_col_width            = 30;
      const unsigned int ncalls_col_width     = 11;
      const unsigned int tot_time_col_width   = 12;
      const unsigned int avg_time_col_width   = 12;
      const unsigned int tot_time_incl_sub_col_width   = 12;
      const unsigned int avg_time_incl_sub_col_width   = 12;
      const unsigned int pct_active_col_width = 9;
      const unsigned int pct_active_incl_sub_col_width = 9;

      // Iterator to be used to loop over the map of timed events
      std::map<std::pair<std::string,std::string>, PerfData>::const_iterator pos;

      // Reset the event column width based on the longest event name plus
      // a possible 2-character indentation, plus a space.
      for (pos = log.begin(); pos != log.end(); ++pos)
        if (pos->first.second.size()+3 > event_col_width)
          event_col_width = cast_int<unsigned int>
            (pos->first.second.size()+3);

      // Set the total width of the column
      const unsigned int total_col_width =
        event_col_width     +
        ncalls_col_width    +
        tot_time_col_width  +
        avg_time_col_width  +
        tot_time_incl_sub_col_width  +
        avg_time_incl_sub_col_width  +
        pct_active_col_width+
        pct_active_incl_sub_col_width+1;

      // Print dashed line
      oss << ' '
          << std::string(total_col_width, '-')
          << '\n';

      {
        // Construct temporary message string
        std::ostringstream temp;
        temp << "| " << label_name << " Performance: Alive time=" << elapsed_time
             << ", Active time=" << total_time;

        // Get the size of the temporary string
        const unsigned int temp_size = cast_int<unsigned int>
          (temp.str().size());

        // Send the temporary message to the output
        oss << temp.str();

        // If this string is longer than the previously computed total
        // column width, skip the additional formatting... this shouldn't
        // happen often, hopefully.  Add two additional characters for a
        // space and a "|" character at the end.
        if (temp_size < total_col_width+2)
          oss << std::setw(total_col_width - temp_size + 2)
              << std::right
              << "|";

        oss << '\n';
      }

      // Print dashed line
      oss << ' '
          << std::string(total_col_width, '-')
          << '\n';


      // Write out the header for the events listing
      oss << "| "
          << std::setw(event_col_width)
          << std::left
          << "Event"
          << std::setw(ncalls_col_width)
          << std::left
          << "nCalls"
          << std::setw(tot_time_col_width)
          << std::left
          << "Total Time"
          << std::setw(avg_time_col_width)
          << std::left
          << "Avg Time"
          << std::setw(tot_time_incl_sub_col_width)
          << std::left
          << "Total Time"
          << std::setw(avg_time_incl_sub_col_width)
          << std::left
          << "Avg Time"
          << std::setw(pct_active_col_width+pct_active_incl_sub_col_width)
          << std::left
          << "% of Active Time"
          << "|\n"
          << "| "
          << std::setw(event_col_width)
          << std::left
          << ""
          << std::setw(ncalls_col_width)
          << std::left
          << ""
          << std::setw(tot_time_col_width)
          << std::left
          << "w/o Sub"
          << std::setw(avg_time_col_width)
          << std::left
          << "w/o Sub"
          << std::setw(tot_time_incl_sub_col_width)
          << std::left
          << "With Sub"
          << std::setw(avg_time_incl_sub_col_width)
          << std::left
          << "With Sub"
          << std::setw(pct_active_col_width)
          << std::left
          << "w/o S"
          << std::setw(pct_active_incl_sub_col_width)
          << std::left
          << "With S"
          << "|\n|"
          << std::string(total_col_width, '-')
          << "|\n|"
          << std::string(total_col_width, ' ')
          << "|\n";

      unsigned int summed_function_calls = 0;
      double       summed_total_time     = 0;
      double       summed_percentage     = 0;

      std::string last_header("");

      for (pos = log.begin(); pos != log.end(); ++pos)
        {
          const PerfData & perf_data = pos->second;

          // Only print the event if the count is non-zero.
          if (perf_data.count != 0)
            {
              const unsigned int perf_count    = perf_data.count;
              const double       perf_time     = perf_data.tot_time;
              const double       perf_avg_time = perf_time / static_cast<double>(perf_count);
              const double       perf_time_incl_sub     = perf_data.tot_time_incl_sub;
              const double       perf_avg_time_incl_sub = perf_time_incl_sub / static_cast<double>(perf_count);
              const double       perf_percent  = (total_time != 0.) ? perf_time / total_time * 100. : 0.;
              const double       perf_percent_incl_sub  = (total_time != 0.) ? perf_time_incl_sub / total_time * 100. : 0.;

              summed_function_calls += perf_count;
              summed_total_time     += perf_time;
              summed_percentage     += perf_percent;

              // Print the event name
              if (pos->first.first == "")
                oss << "| "
                    << std::setw(event_col_width)
                    << std::left
                    << pos->first.second;

              else
                {
                  if (last_header != pos->first.first)
                    {
                      last_header = pos->first.first;

                      // print blank line followed by header name
                      // (account for additional space before the
                      // header)
                      oss << "|"
                          << std::string(total_col_width, ' ')
                          << "|\n| "
                          << std::setw(total_col_width-1)
                          << std::left
                          << pos->first.first
                          << "|\n";
                    }

                  oss << "|   "
                      << std::setw(event_col_width-2)
                      << std::left
                      << pos->first.second;
                }


              // Print the number of calls to the event.
              oss << std::setw(ncalls_col_width)
                  << perf_count;

              // Save the original stream flags
              std::ios_base::fmtflags out_flags = oss.flags();

              // Print the total time spent in the event
              oss << std::fixed
                  << std::setprecision(4)
                  << std::setw(tot_time_col_width)
                  << std::left
                  << perf_time;


              // Print the average time per function call
              oss << std::fixed
                  << std::setprecision(6)
                  << std::setw(avg_time_col_width)
                  << std::left
                  << perf_avg_time;

              // Print the total time spent in the event incl. sub-events
              oss << std::fixed
                  << std::setprecision(4)
                  << std::setw(tot_time_incl_sub_col_width)
                  << std::left
                  << perf_time_incl_sub;

              // Print the average time per function call incl. sub-events
              oss << std::fixed
                  << std::setprecision(6)
                  << std::setw(avg_time_incl_sub_col_width)
                  << std::left
                  << perf_avg_time_incl_sub;

              // Print the percentage of the time spent in the event
              oss << std::fixed
                  << std::setprecision(2)
                  << std::setw(pct_active_col_width)
                  << std::left
                  << perf_percent;

              // Print the percentage of the time spent in the event incl. sub-events
              oss << std::fixed
                  << std::setprecision(2)
                  << std::setw(pct_active_incl_sub_col_width)
                  << std::left
                  << perf_percent_incl_sub;

              // Reset the stream flags
              oss.flags(out_flags);

              oss << "|\n";
            }
        }

      oss << ' '
          << std::string(total_col_width, '-')
          << "\n| "
          << std::setw(event_col_width)
          << std::left
          << "Totals:";

      // Print the total number of logged function calls
      // For routines which are called many times, summed_function_calls may
      // exceed 7 digits.  If this happens use, scientific notation.
      if (summed_function_calls < 9999999)
        oss << std::setw(ncalls_col_width)
            << summed_function_calls;

      else
        {
          // Save the original stream flags
          std::ios_base::fmtflags out_flags = oss.flags();

          oss << std::scientific
              << std::setprecision(3)
              << std::setw(ncalls_col_width)
              << std::left
              << static_cast<Real>(summed_function_calls);

          // Reset the stream flags
          oss.flags(out_flags);
        }

      // Print the total time spent in logged function calls.  Don't bother saving/restoring
      // the flags here since we are almost done with this stream anyway...
      oss << std::fixed
          << std::setprecision(4)
          << std::setw(tot_time_col_width)
          << std::left
          << summed_total_time;

      // Null, the average time doesn't make sense as a total
      oss << std::setw(avg_time_col_width) << "";

      // Same for times that include sub-events
      oss << std::setw(tot_time_incl_sub_col_width)
          << ""
          << std::setw(avg_time_incl_sub_col_width)
          << "";

      // Print the total percentage followed by dashed line
      oss << std::fixed
          << std::setprecision(2)
          << std::setw(pct_active_col_width)
          << std::left
          << summed_percentage
          << std::setw(pct_active_incl_sub_col_width)
          << ""
          << "|\n "
          << std::string(total_col_width, '-')
          << '\n';
    }

  return oss.str();
}



std::string PerfLog::get_log() const
{
  std::ostringstream oss;

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
              oss << get_info_header();
            }
          oss << get_perf_info();
        }
    }

  return oss.str();
}



void PerfLog::print_log() const
{
  if (log_events)
    {
      // Check to see if the log_string is empty, and if so,
      // avoid printing an unnecessary newline.
      std::string log_string = this->get_log();
      if (log_string.size() > 0)
        libMesh::out << log_string << std::endl;
    }
}

PerfData PerfLog::get_perf_data(const std::string & label, const std::string & header)
{
  return log[std::make_pair(header, label)];
}

void PerfLog::start_event(const std::string & label,
                          const std::string & header)
{
  this->push(label,header);
}



void PerfLog::stop_event(const std::string & label,
                         const std::string & header)
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



void PerfLog::split_on_whitespace(const std::string & input, std::vector<std::string> & output) const
{
  // Check for easy return
  if (input.size()==0)
    return;

  // Here we hard-code the string to split on, since the algorithm below
  // is somewhat specific to it...
  const std::string split_on("' '");

  size_t current_pos = 0;
  while (true)
    {
      // Find next end location
      size_t end_pos = input.find(split_on, current_pos);

      if (end_pos != std::string::npos)
        {
          // Create substring.  Note: the second argument to substr is
          // the *length* of string to create, not the ending position!
          output.push_back( input.substr(current_pos, end_pos - current_pos + 1) );

          // Update search starting position, make sure to go past the end of the split_on string, but
          // include the previous single quote (hence the -1).
          current_pos = end_pos + split_on.size() - 1;
        }
      else
        {
          // Push back whatever remains of the string onto the output.
          // Note that substr with only 1 argument pushes back
          // whatever remains of the string.  This also handles the
          // case where the string does not contain any matches.
          output.push_back( input.substr(current_pos) );

          // We are done searching the string, so break out of the while loop
          break;
        }
    }
}

} // namespace libMesh
