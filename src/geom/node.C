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


processor_id_type
Node::choose_processor_id(processor_id_type pid1, processor_id_type pid2) const
{
  if (pid1 == DofObject::invalid_processor_id)
    return pid2;

  // Do we want the new load-balanced node partitioning heuristic
  // instead of the default partitioner-friendlier heuristic?
  static bool load_balanced_nodes =
    libMesh::on_command_line ("--load-balanced-nodes");

  // For better load balancing, we can use the min
  // even-numberered nodes and the max for odd-numbered.
  if (load_balanced_nodes)
    {
      if (this->id() % 2 &&
          pid2 != DofObject::invalid_processor_id)
        return std::max(pid1, pid2);
      else
        return std::min(pid1, pid2);
    }

  // Our default behavior, which puts too many nodes on lower MPI
  // ranks but which keeps elements' nodes on the same partition more
  // often, is simply:
  return std::min(pid1, pid2);
}


} // namespace libMesh
