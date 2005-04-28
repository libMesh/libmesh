// $Id: mesh_triangle_support.h,v 1.1 2005-04-28 20:38:21 jwpeterson Exp $
 
// The libMesh Finite Element Library.
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


#ifndef __mesh_triangle_support_h__
#define __mesh_triangle_support_h__

#include "libmesh_config.h"

#ifdef HAVE_TRIANGLE

// Note: libmesh_common.h defines REAL, which is required by triangle.
// Therefore, we need to include it first.
#include "libmesh_common.h"

extern "C" {
#include "triangle.h"
}







#endif // HAVE_TRIANGLE

#endif
