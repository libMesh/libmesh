// $Id: node.C,v 1.11 2007-10-21 20:48:45 benkirk Exp $

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




// Local includes
#include "node.h"




// ------------------------------------------------------------
// Node class static member initialization
//const unsigned int Node::invalid_id = libMesh::invalid_uint;


bool Node::operator==(const DofObject& rhs) const
{
  // Cast rhs to a Node*
  const Node* rhs_node = dynamic_cast<const Node*>(&rhs);

  // If we can't cast to a Node* then rhs must be an Elem
  if(rhs_node == NULL)
    return false;

  // Explicitly calling the operator== defined in Point
  return this->Point::operator==(*rhs_node);
}
