// $Id: tree.C,v 1.6 2003-02-13 22:56:14 benkirk Exp $

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

// Local includes
#include "tree.h"



// ------------------------------------------------------------
// Tree class method
template <unsigned int N>
Elem* Tree<N>::find_element(const Point& p) const
{
  return root.find_element(p);
}


// ------------------------------------------------------------
// Explicit Instantiations
template class Tree<4>;
template class Tree<8>;









