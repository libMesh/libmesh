// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <sstream>

// Local includes
#include "libmesh/reference_counter.h"

namespace libMesh
{



// ------------------------------------------------------------
// ReferenceCounter class static member initializations
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

ReferenceCounter::Counts ReferenceCounter::_counts;

#endif

bool ReferenceCounter::_enable_print_counter = true;
Threads::atomic<unsigned int> ReferenceCounter::_n_objects;
Threads::spin_mutex  ReferenceCounter::_mutex;


// ------------------------------------------------------------
// ReferenceCounter class members
std::string ReferenceCounter::get_info ()
{
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

  std::ostringstream oss;

  oss << '\n'
      << " ---------------------------------------------------------------------------- \n"
      << "| Reference count information                                                |\n"
      << " ---------------------------------------------------------------------------- \n";

  for (const auto & pr : _counts)
    {
      const std::string name(pr.first);
      const unsigned int creations    = pr.second.first;
      const unsigned int destructions = pr.second.second;

      oss << "| " << name << " reference count information:\n"
          << "|  Creations:    " << creations    << '\n'
          << "|  Destructions: " << destructions << '\n';
    }

  oss << " ---------------------------------------------------------------------------- \n";

  return oss.str();

#else

  return "";

#endif
}





// avoid unused variable warnings
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)

void ReferenceCounter::print_info (std::ostream & out_stream)
{
  if (_enable_print_counter)
    out_stream << ReferenceCounter::get_info();
}

#else

void ReferenceCounter::print_info (std::ostream & /* out_stream */)
{}

#endif

void ReferenceCounter::enable_print_counter_info()
{
  _enable_print_counter = true;
  return;
}

void ReferenceCounter::disable_print_counter_info()
{
  _enable_print_counter = false;
  return;
}

} // namespace libMesh
