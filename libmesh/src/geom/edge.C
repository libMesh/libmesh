// $Id$

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
#include "edge.h"
#include "node_elem.h"


AutoPtr<DofObject> Edge::side (const unsigned int i) const
{
  assert(i < 2);
  const Elem* parent = this;
  Elem *nodeelem = new NodeElem(const_cast<Elem*>(parent));
  nodeelem->set_node(0) = this->get_node(i);
  return AutoPtr<DofObject>(nodeelem);
}


AutoPtr<Elem> Edge::build_side (const unsigned int i, bool) const
{
  assert(i < 2);
  const Elem* parent = this;
  Elem *nodeelem = new NodeElem(const_cast<Elem*>(parent));
  nodeelem->set_node(0) = this->get_node(i);
  return AutoPtr<Elem>(nodeelem);
}
