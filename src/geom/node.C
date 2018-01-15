// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/node.h"

namespace libMesh
{




// ------------------------------------------------------------
// Node class static member initialization
//const unsigned int Node::invalid_id = libMesh::invalid_uint;


bool Node::operator==(const Node & rhs) const
{
  // Explicitly calling the operator== defined in Point
  return this->Point::operator==(rhs);
}



void Node::print_info (std::ostream & os) const
{
  os << this->get_info()
     << std::endl;
}



std::string Node::get_info () const
{
  std::ostringstream oss;

  oss << "  Node id()=";

  if (this->valid_id())
    oss << this->id();
  else
    oss << "invalid";

  oss << ", processor_id()=" << this->processor_id() <<
    ", Point=" << *static_cast<const Point *>(this) << '\n';

  oss << "    DoFs=";
  for (unsigned int s=0; s != this->n_systems(); ++s)
    for (unsigned int v=0; v != this->n_vars(s); ++v)
      for (unsigned int c=0; c != this->n_comp(s,v); ++c)
        oss << '(' << s << '/' << v << '/' << this->dof_number(s,v,c) << ") ";

  return oss.str();
}


} // namespace libMesh
