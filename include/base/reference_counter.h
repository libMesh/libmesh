// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_REFERENCE_COUNTER_H
#define LIBMESH_REFERENCE_COUNTER_H

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/threads.h"
#include "libmesh/libmesh.h" // libMesh::on_command_line

// C++ includes
#include <iostream>
#include <string>
#include <map>

namespace libMesh
{



/**
 * This is the base class for enabling reference counting.  It
 * should not be used by the user, thus it has a private constructor.
 *
 * \author Benjamin S. Kirk, 2002-2007
 */

// ------------------------------------------------------------
// ReferenceCounter class definition
class ReferenceCounter
{
protected:

  /**
   * Constructor. Protected so that you cannont
   * instantiate a \p ReferenceCounter, only derive
   * from it.
   */
  ReferenceCounter ();

public:

  /**
   * Destructor.
   */
  ~ReferenceCounter ();

  /**
   * Gets a string containing the reference information.
   */
  static std::string get_info ();

  /**
   * Prints the reference information, by default to \p libMesh::out.
   */
  static void print_info (std::ostream& out = libMesh::out);

  /**
   * Prints the number of outstanding (created, but not yet
   * destroyed) objects.
   */
  static unsigned int n_objects ()
  { return _n_objects; }

  /**
   * Methods to enable/disable the reference counter output
   * from print_info()
   */
  static void enable_print_counter_info();
  static void disable_print_counter_info();


protected:

#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  /**
   * Increments the construction counter. Should be called in
   * the constructor of any derived class that will be
   * reference counted.
   */
  void increment_constructor_count (const std::string& name);

  /**
   * Increments the destruction counter. Should be called in
   * the destructor of any derived class that will be
   * reference counted.
   */
  void increment_destructor_count (const std::string& name);

  /**
   * Data structure to log the information.  The log is
   * identified by the class name.
   */
  typedef std::map<std::string, std::pair<unsigned int,
                                          unsigned int> > Counts;

  /**
   * Actually holds the data.
   */
  static Counts _counts;

#endif

  /**
   * The number of objects.  Print the reference count
   * information when the number returns to 0.
   */
  static Threads::atomic<unsigned int> _n_objects;

  /**
   * Mutual exclusion object to enable thread-safe reference counting.
   */
  static Threads::spin_mutex _mutex;

  /**
   * Flag to control whether reference count information
   * is printed when print_info is called.
   */
  static bool _enable_print_counter;
};



// ------------------------------------------------------------
// ReferenceCounter class inline methods
inline ReferenceCounter::ReferenceCounter()
{
  ++_n_objects;
}



inline ReferenceCounter::~ReferenceCounter()
{
  --_n_objects;
}






#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
inline
void ReferenceCounter::increment_constructor_count (const std::string& name)
{
  Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
  std::pair<unsigned int, unsigned int>& p = _counts[name];

  p.first++;
}
#endif



#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
inline
void ReferenceCounter::increment_destructor_count (const std::string& name)
{
  Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
  std::pair<unsigned int, unsigned int>& p = _counts[name];

  p.second++;
}
#endif


} // namespace libMesh


#endif // LIBMESH_REFERENCE_COUNTER_H
