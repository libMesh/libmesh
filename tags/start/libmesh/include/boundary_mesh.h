// $Id: boundary_mesh.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __boundary_mesh_h__
#define __boundary_mesh_h__

// C++ Includes   -----------------------------------
#include <iostream>
#include <vector>
#include <string>



// Local Includes -----------------------------------
#include "mesh_base.h"



/**
 * The \p BoundaryMesh is a \p Mesh in its own right, but it
 * contains a description of the boundary of some other mesh.
 * This is useful for writing the boundary of a domain for inspecting
 * boundary conditions and other things.
 */

// ------------------------------------------------------------
// BoundaryMesh class definition
class BoundaryMesh : public MeshBase
{
 public:
  
  BoundaryMesh(unsigned int d,
	       unsigned int proc_id=0);
  
  ~BoundaryMesh();

 private:
};




#endif
