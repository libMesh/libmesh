// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


// Local includes
#include "libmesh_common.h"
#include "o_string_stream.h"

// C++ includes
#include <string>
#include <map>
#include <sys/time.h>

#ifdef HAVE_LOCALE
#include <locale>
#endif

// Forward Declarations
// class OStringStream;


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
    {}

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
   * Constructor.  \p label_name is the name of the object, which
   * will bw print in the log to distinguish it from other objects.
   * \p log_events is a flag to optionally
   * disable logging.  You can use this flag to turn off
   * logging without touching any other code.
   */
  PerfLog(const std::string& label_name="",
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
  void start_event(const std::string &label,
		   const std::string &header="");

  /**
   * Stop monitoring the event named \p label.
   */
  void stop_event(const std::string &label,
		  const std::string &header="");

  /**
   * Suspend monitoring of the event. 
   */
  void pause_event(const std::string &label,
		   const std::string &header="");

  /**
   * Restart monitoring the event.
   */
  void restart_event(const std::string &label,
		     const std::string &header="");
  
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

  /**
   * @returns the total time spent on this event.
   */
  double get_total_time() const
    {return total_time;}
 
   
 private:

  
  /**
   * The label for this object.
   */
  const std::string label_name;

  /**
   * Flag to optionally disable all logging.
   */
  const bool log_events;

  /**
   * The total running time for recorded events.
   */  
  double total_time;
  
  /**
   * The time we were constructed or last cleared.
   */
  struct timeval tstart;

  /**
   * The actual log.
   */
  std::map<std::pair<std::string,
		     std::string>,
	   PerfData> log;
  
  /**
   * Flag indicating if print_log() has been called.
   * This is used to print a header with machine-specific
   * data the first time that print_log() is called.
   */
  static bool called;

  /**
   * Prints a line of 'n' repeated characters 'c'
   * to the output string stream "out".
   */
  void _character_line(const unsigned int n,
		       const char c,
		       OStringStream& out) const;
};



// ------------------------------------------------------------
// PerfLog class inline member funcions

// Typedefs we might need
#ifdef HAVE_LOCALE
typedef std::ostreambuf_iterator<char, std::char_traits<char> > TimeIter;
typedef std::time_put<char, TimeIter> TimePut;
#endif

#endif
