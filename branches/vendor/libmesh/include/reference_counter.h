 // $Id: reference_counter.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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

// C++ includes
#include <iostream>
#include <string>
#include <map>

// Local includes
#include "mesh_config.h"
#include "mesh_common.h"




/**
 * This class implements reference counting. Any class that
 * is derived from this class that properly implements the
 * \p class_name() member will get reference counted, provided
 * that the library is configured with \p --enable-reference-counting.
 * The derived class must call \p increment_constructor_count()
 * in its constructor and \p increment_destructor_count()
 * in its destructor.
 * \par
 * If the library is configured with \p --disable-reference-counting
 * then this class does nothing.  All members are inlined and empty,
 * so they should effectively disappear.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.1.1.1 $
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
  ReferenceCounter () {};

public:
  
  /**
   * Destructor.
   */
  virtual ~ReferenceCounter () {};

  /**
   * @returns the class name as a \p std::string.  This method
   * must be implemented in the derived class.
   */
  virtual std::string class_name () const = 0;

  /**
   * Increments the construction counter. Should be called in
   * the constructor of any derived class that will be
   * reference counted.
   */
  void increment_constructor_count ();
  
  /**
   * Increments the destruction counter. Should be called in
   * the destructor of any derived class that will be
   * reference counted.
   */ 
  void increment_destructor_count  ();
  
  /**
   * Prints the reference information to \p std::cout.
   */
  static void print_info ();

  
private:

  
#ifdef ENABLE_REFERENCE_COUNTING
  
  /**
   * Data structure to log the information.  The log is
   * identified by the class name.
   */
  typedef std::map<std::string, std::pair<unsigned int,
					  unsigned int> > Counts;

  /**
   * Actually holds the data.
   */
  static Counts counts;
  
#endif
  
};



// ------------------------------------------------------------
// ReferenceCounter class inline methods
inline
void ReferenceCounter::increment_constructor_count ()
{
#ifdef ENABLE_REFERENCE_COUNTING

  std::pair<unsigned int, unsigned int>& p = counts[class_name()];

  p.first++;

#endif
};



inline
void ReferenceCounter::increment_destructor_count ()
{
#ifdef ENABLE_REFERENCE_COUNTING

  std::pair<unsigned int, unsigned int>& p = counts[class_name()];

  p.second++;

#endif
};



inline
void ReferenceCounter::print_info ()
{
#ifdef ENABLE_REFERENCE_COUNTING
  
  for (Counts::iterator it = counts.begin();
       it != counts.end(); ++it)
    {
      const std::string& name         = it->first;
      const unsigned int creations    = it->second.first;
      const unsigned int destructions = it->second.second;

      std::cout << name
		<< " class reference count information:"
		<< std::endl
		<< " Creations:    " << creations
		<< std::endl
		<< " Destructions: " << destructions
		<< std::endl;

      if (creations != destructions)
	std::cout << " WARNING: class "
		  << name << " "
		  << creations - destructions
		  << " items leaked!"
		  << std::endl;
      
    };

#endif
};

#endif




