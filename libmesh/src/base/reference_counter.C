// $Id: reference_counter.C,v 1.5 2003-02-10 03:55:51 benkirk Exp $

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
#include <sstream>

// Local includes
#include "reference_counter.h"



// ------------------------------------------------------------
// ReferenceCounter class static member initializations
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

ReferenceCounter::Counts ReferenceCounter::_counts;
unsigned int             ReferenceCounter::_n_objects=0;

#endif



// ------------------------------------------------------------
// ReferenceCounter class members
std::string ReferenceCounter::get_info ()
{
  std::ostringstream out;
  
#if defined(ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  out << std::endl
      << " ---------------------------------------------------------------------- "  << std::endl
      << "| Reference count information                                          |" << std::endl
      << " ---------------------------------------------------------------------- "  << std::endl;
  
  for (Counts::iterator it = _counts.begin();
       it != _counts.end(); ++it)
    {
      const std::string name(it->first);
      const unsigned int creations    = it->second.first;
      const unsigned int destructions = it->second.second;

      out << "| "
	  << name
	  << " class reference count information:"
	  << std::endl
	  << "| Creations:    " << creations
	  << std::endl
	  << "| Destructions: " << destructions
	  << std::endl;

      if (creations != destructions)
	out << "| WARNING: class "
	    << name << " "
	    << creations - destructions
	    << " items leaked!"
	    << std::endl;
      
    };
  
  out << " ---------------------------------------------------------------------- "  << std::endl;

#endif

  return out.str();
};
