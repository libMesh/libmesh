// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __elem_range_h__
#define __elem_range_h__

// C++ includes

// Local includes
#include "mesh_base.h"
#include "elem.h"
#include "stored_range.h"

namespace libMesh
{

typedef StoredRange<MeshBase::element_iterator,             Elem*>      ElemRange;
typedef StoredRange<MeshBase::const_element_iterator, const Elem*> ConstElemRange;

} // namespace libMesh

#endif // end #ifndef __elem_range_h__
