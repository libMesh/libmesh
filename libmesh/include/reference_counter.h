 // $Id: reference_counter.h,v 1.2 2003-01-20 16:31:29 jwpeterson Exp $

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



#ifndef __reference_counter_h__
#define __reference_counter_h__

// Local includes
#include "mesh_config.h"


// C++ includes
#ifdef ENABLE_REFERENCE_COUNTING
#include <iostream>
#include <string>
#include <map>
#endif



/**
 * This is the base class for enabling reference counting.  It
 * should not be used by the user, thus it has a private constructor.
 *
 * \author Benjamin S. Kirk, 2002-2003
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
  virtual ~ReferenceCounter ();

  /**
   * Gets a string containing the reference information.
   */
  static std::string get_info ();
  
  /**
   * Prints the reference information to \p std::cout.
   */
  static void print_info ();

  
protected:

  
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


#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
  
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

  /**
   * The number of objects.  Print the reference count
   * information when the number returns to 0.
   */
  static unsigned int _n_objects;
#endif
};



// ------------------------------------------------------------
// ReferenceCounter class inline methods
inline ReferenceCounter::ReferenceCounter()
{
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
  
  _n_objects++;

#endif
};



inline ReferenceCounter::~ReferenceCounter()
{
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
  
  _n_objects--;

  if (_n_objects == 0)
    print_info();

#endif
};



inline
void ReferenceCounter::print_info ()
{
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  std::cout << get_info();

#endif
};



inline
void ReferenceCounter::increment_constructor_count (const std::string& name)
{
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  std::pair<unsigned int, unsigned int>& p = _counts[name];

  p.first++;

#endif
};



inline
void ReferenceCounter::increment_destructor_count (const std::string& name)
{
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  std::pair<unsigned int, unsigned int>& p = _counts[name];

  p.second++;

#endif
};



#endif // end #ifndef __reference_counter_h__




