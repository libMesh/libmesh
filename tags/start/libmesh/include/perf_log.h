// $Id: perf_log.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __perflog_h__
#define __perflog_h__

#include "mesh_config.h"


// C++ includes
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <sys/time.h>

#ifdef HAVE_LOCALE
#include <locale>
#endif

// Local includes
#include "mesh_common.h"

/**
 * The \p PerfData class simply contains the performance
 * data that is recorded for individual events.
 */

// ------------------------------------------------------------
// PerfData class definition
class PerfData
{
 public:

  /**
   * Constructor.  Initializes data to be empty.
   */
  PerfData () :
    tot_time(0.),
    count(0),
    open(false)
    {};

  /**
   * Total time spent in this event.
   */
  double tot_time;

  /**
   * Structure defining when the event
   * was last started.
   */
  struct timeval tstart;

  /**
   * The number of times this event has
   * been executed
   */
  unsigned int count;

  /**
   * Flag indicating if we are currently
   * monitoring this event.  Should only
   * be true while the event is executing.
   */
  bool open;
};




/**
 * The \p PerfLog class allows monitoring of specific events.
 * An event is defined by a unique string that functions as
 * a label.  Each time the event is executed data are recorded.
 * This class is particulary useful for finding performance
 * bottlenecks. 
 *
 */

// ------------------------------------------------------------
// PerfLog class definition
class PerfLog
{

 public:

  /**
   * Constructor.  \p class_name is the name of the class
   * we are monitoring.  \p log_events is a flag to optionally
   * disable logging.  You can use this flag to turn off
   * logging without touching any other code.
   */
  PerfLog(std::string class_name="",
	  const bool log_events=true);

  /**
   * Destructor. Calls \p clear() and \p print_log().
   */
  ~PerfLog();
  
  /**
   * Clears all the internal data and returns the
   * data structures to a pristine state.  This function
   * checks to see if it is currently monitoring any
   * events, and if so errors.  Be sure you are not
   * logging any events when you call this function.
   */
  void clear();

  /**
   * Start monitoring the event named \p label.
   */
  void start_event(const std::string &label);

  /**
   * Stop monitoring the event named \p label.
   */
  void stop_event(const std::string &label);

  /**
   * Suspend monitoring of the event. 
   */
  void pause_event(const std::string &label);

  /**
   * Restart monitoring the event.
   */
  void restart_event(const std::string &label);
  
  /**
   * @returns a string containing:
   * (1) Basic machine information (if first call)
   * (2) The performance log
   */
  std::string get_log() const;
  
  /**
   * @returns a string containing ONLY the information header.
   */
  std::string get_info_header() const;

  /**
   * @returns a string containing ONLY the log information
   */
  std::string get_perf_info() const;
  
  /**
   * Print the log.
   */
  void print_log() const;

    
 private:

  
  /**
   * The name of the class we are monitoring.
   */
  const std::string class_name;

  /**
   * Flag to optionally disable all logging.
   */
  const bool log_events;

  /**
   * The time we were constructed or last cleared.
   */
  struct timeval tstart;

  /**
   * The actual log.
   */
  std::map<std::string, PerfData> log;

  /**
   * Flag indicating if print_log() has been called.
   * This is used to print a header with machine-specific
   * data the first time that print_log() is called.
   */
  static bool called;
};



// ------------------------------------------------------------
// PerfLog class inline member funcions
inline
void PerfLog::start_event(const std::string &label)
{
  if (log_events)
    {
      // make sure we aren't currently
      // monitoring this event
      if (log[label].open)
	{
	  std::cerr << "ERROR logging event " << label << std::endl
		    << "Did you forget to stop logging it?" << std::endl;

	  error();
	}

      log[label].open = true;
      
      gettimeofday (&log[label].tstart, NULL);
    }
};



inline
void PerfLog::stop_event(const std::string &label)
{
  if (log_events)
    {
      // make sure we are currently
      // monitoring this event
      if (!log[label].open)
	{
	  std::cerr << "ERROR logging event " << label << std::endl
		    << "Did you forget to start logging it?" << std::endl;

	  error();
	}
      
      log[label].open = false;
      
      struct timeval tstop;

      gettimeofday (&tstop, NULL);

      const double elapsed_time = ((double) (tstop.tv_sec  - log[label].tstart.tv_sec)) +
	                          ((double) (tstop.tv_usec - log[label].tstart.tv_usec))/1000000.;

      log[label].tot_time += elapsed_time;
      log[label].count++;	 
    }
};



inline
void PerfLog::pause_event(const std::string &label)
{
  if (log_events)
    {
      // make sure we are currently
      // monitoring this event
      if (!log[label].open)
	{
	  std::cerr << "ERROR pausing event " << label << std::endl
		    << "Did you forget to start logging it?" << std::endl;
	  
	  error();
	}
      
      struct timeval tstop;

      gettimeofday (&tstop, NULL);

      const double elapsed_time = ((double) (tstop.tv_sec  - log[label].tstart.tv_sec)) +
	                          ((double) (tstop.tv_usec - log[label].tstart.tv_usec))/1000000.;

      log[label].tot_time += elapsed_time;
    }
};



inline
void PerfLog::restart_event(const std::string &label)
{
  if (log_events)
    {
      // make sure we are currently
      // monitoring this event
      if (!log[label].open)
	{
	  std::cerr << "ERROR restarting event " << label << std::endl
		    << "Did you forget to start logging it?" << std::endl;
	  
	  error();
	}
      
      gettimeofday (&log[label].tstart, NULL);
    }
};

// Typedefs we might need
#ifdef HAVE_LOCALE
typedef std::ostreambuf_iterator<char, std::char_traits<char> > TimeIter;
typedef std::time_put<char, TimeIter> TimePut;
#endif

#endif
