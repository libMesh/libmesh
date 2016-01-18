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



#ifndef LIBMESH_TIMESTAMP_H
#define LIBMESH_TIMESTAMP_H

// C++ includes
#include <iostream>
#include <string>

namespace libMesh
{

namespace Utility
{
// Two helpful utility programs for getting and printing a timestamp in
// the format:
// 04 February 2010 01:56:43 PM
std::string get_timestamp();

inline void print_timestamp (std::ostream & target = std::cout)
{
  target << get_timestamp() << std::endl;
}

}

} // namespace libMesh


#endif // LIBMESH_TIMESTAMP_H
