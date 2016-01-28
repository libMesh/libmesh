// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PERFLOG_H
#define LIBMESH_PERFLOG_H


// Local includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <cstddef>
#include <map>
#include <stack>
#include <string>
#include <vector>
#include <sys/time.h>

namespace libMesh
{

/**
 * The \p PerfData class simply contains the performance
 * data that is recorded for individual events.
 *
 * \author Benjamin Kirk
 * \date 2003
 */
class PerfData
{
public:

  /**
   * Constructor.  Initializes data to be empty.
   */
  PerfData () :
    tot_time(0.),
    tot_time_incl_sub(0.),
    tstart(),
    tstart_incl_sub(),
    count(0),
    open(false),
    called_recursively(0)
  {}


  /**
   * Total time spent in this event.
   */
  double tot_time;

  /**
   * Total time spent in this event, including sub-events.
   */
  double tot_time_incl_sub;

  /**
   * Structure defining when the event
   * was last started.
   */
  struct timeval tstart;

  /**
   * Structure defining when the event
   * was last started, including sub-events.
   */
  struct timeval tstart_incl_sub;

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

  void   start ();
  void   restart ();
  double pause ();
  double stopit ();

  int called_recursively;

protected:
  double stop_or_pause(const bool do_stop);
};




/**
 * The \p PerfLog class allows monitoring of specific events.
 * An event is defined by a unique string that functions as
 * a label.  Each time the event is executed data are recorded.
 * This class is particulary useful for finding performance
 * bottlenecks.
 *
 * \author Benjamin Kirk
 * \date 2003
 */
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
  PerfLog(const std::string & label_name="",
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
   * Disables performance logging for an active object.
   */
  void disable_logging() { log_events = false; }

  /**
   * Enables performance logging for an active object.
   */
  void enable_logging() { log_events = true; }

  /**
   * Returns true iff performance logging is enabled
   */
  bool logging_enabled() const { return log_events; }

  /**
   * Push the event \p label onto the stack, pausing any active event.
   */
  void push (const char * label,
             const char * header="");

  /**
   * Push the event \p label onto the stack, pausing any active event.
   */
  void push (const std::string & label,
             const std::string & header="");

  /**
   * Pop the event \p label off the stack, resuming any lower event.
   */
  void pop (const char * label,
            const char * header="");

  /**
   * Pop the event \p label off the stack, resuming any lower event.
   */
  void pop (const std::string & label,
            const std::string & header="");

  /**
   * Start monitoring the event named \p label.
   */
  void start_event(const std::string & label,
                   const std::string & header="");

  /**
   * Stop monitoring the event named \p label.
   */
  void stop_event(const std::string & label,
                  const std::string & header="");

  /**
   * Suspend monitoring of the event.
   */
  void pause_event(const std::string & label,
                   const std::string & header="");

  /**
   * Restart monitoring the event.
   */
  void restart_event(const std::string & label,
                     const std::string & header="");

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
  double get_elapsed_time() const;

  /**
   * @returns the active time
   */
  double get_active_time() const;

  /**
   * Return the PerfData object associated with a label and header.
   */
  PerfData get_perf_data(const std::string & label, const std::string & header="");

private:


  /**
   * The label for this object.
   */
  const std::string label_name;

  /**
   * Flag to optionally disable all logging.
   */
  bool log_events;

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
   * A stack to hold the current performance log trace.
   */
  std::stack<PerfData*> log_stack;

  /**
   * Flag indicating if print_log() has been called.
   * This is used to print a header with machine-specific
   * data the first time that print_log() is called.
   */
  static bool called;

  /**
   * Splits a string on whitespace into a vector of separate strings.  This is used to make the
   * LIBMESH_CONFIGURE_INFO a little more manageable.
   */
  void split_on_whitespace(const std::string & input,
                           std::vector<std::string> & output) const;
};



// ------------------------------------------------------------
// PerfData class member funcions
inline
void PerfData::start ()
{
  this->count++;
  this->called_recursively++;
  gettimeofday (&(this->tstart), libmesh_nullptr);
  this->tstart_incl_sub = this->tstart;
}



inline
void PerfData::restart ()
{
  gettimeofday (&(this->tstart), libmesh_nullptr);
}



inline
double PerfData::pause ()
{
  return this->stop_or_pause(false);
}


inline
double PerfData::stop_or_pause(const bool do_stop)
{
  // save the start times, reuse the structure we have rather than create
  // a new one.
  const time_t
    tstart_tv_sec  = this->tstart.tv_sec,
    tstart_tv_usec = this->tstart.tv_usec;

  gettimeofday (&(this->tstart), libmesh_nullptr);

  const double elapsed_time = (static_cast<double>(this->tstart.tv_sec  - tstart_tv_sec) +
                               static_cast<double>(this->tstart.tv_usec - tstart_tv_usec)*1.e-6);

  this->tot_time += elapsed_time;

  if(do_stop)
    {
      const double elapsed_time_incl_sub = (static_cast<double>(this->tstart.tv_sec  - this->tstart_incl_sub.tv_sec) +
                                            static_cast<double>(this->tstart.tv_usec - this->tstart_incl_sub.tv_usec)*1.e-6);

      this->tot_time_incl_sub += elapsed_time_incl_sub;
    }

  return elapsed_time;
}



inline
double PerfData::stopit ()
{
  // stopit is just similar to pause except that it decrements the
  // recursive call counter

  this->called_recursively--;
  return this->stop_or_pause(true);
}



// ------------------------------------------------------------
// PerfLog class inline member funcions
inline
void PerfLog::push (const std::string & label,
                    const std::string & header)
{
  if (this->log_events)
    {
      // Get a reference to the event data to avoid
      // repeated map lookups
      PerfData * perf_data = &(log[std::make_pair(header,label)]);

      if (!log_stack.empty())
        total_time +=
          log_stack.top()->pause();

      perf_data->start();
      log_stack.push(perf_data);
    }
}



inline
void PerfLog::push (const char * label,
                    const char * header)
{
  if (this->log_events)
    this->push(std::string(label), std::string(header));
}



inline
void PerfLog::pop (const std::string & libmesh_dbg_var(label),
                   const std::string & libmesh_dbg_var(header))
{
  if (this->log_events)
    {
      libmesh_assert (!log_stack.empty());

#ifndef NDEBUG
      PerfData * perf_data = &(log[std::make_pair(header,label)]);
      if (perf_data != log_stack.top())
        {
          libMesh::err << "PerfLog can't pop (" << header << ',' << label << ')' << std::endl;
          libMesh::err << "From top of stack of running logs:" << std::endl;
          std::map<std::pair<std::string, std::string>, PerfData>::iterator
            i = log.begin(), endi = log.end();
          for (; i != endi; ++i)
            if (&(i->second) == log_stack.top())
              libMesh::err << '(' << i->first.first << ',' << i->first.second << ')' << std::endl;

          libmesh_assert_equal_to (perf_data, log_stack.top());
        }
#endif

      total_time += log_stack.top()->stopit();

      log_stack.pop();

      if (!log_stack.empty())
        log_stack.top()->restart();
    }
}



inline
void PerfLog::pop(const char * label,
                  const char * header)
{
  if (this->log_events)
    this->pop(std::string(label), std::string(header));
}



inline
double PerfLog::get_elapsed_time () const
{
  struct timeval tnow;

  gettimeofday (&tnow, libmesh_nullptr);

  const double elapsed_time = (static_cast<double>(tnow.tv_sec  - tstart.tv_sec) +
                               static_cast<double>(tnow.tv_usec - tstart.tv_usec)*1.e-6);
  return elapsed_time;
}

inline
double PerfLog::get_active_time() const
{
  return total_time;
}

} // namespace libMesh



#endif // LIBMESH_PERFLOG_H
