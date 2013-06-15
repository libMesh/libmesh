// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_ELEM_ASSEMBLY_WITH_CONSTRUCTION_H
#define LIBMESH_ELEM_ASSEMBLY_WITH_CONSTRUCTION_H

// Local includes
#include "libmesh/reference_counted_object.h"

// C++ includes

namespace libMesh
{

/**
 * Extend ElemAssembly to provide access to the RBConstruction
 * class. This is often needed in reduced basis applications,
 * where the RBConstruction object gives us access to the mesh
 * and hence, for example, boundary condition information.
 *
 * @author David J. Knezevic, 2013
 */
class ElemAssemblyWithConstruction : public ElemAssembly
{
public:

  /**
   * Constructor.
   */
  ElemAssemblyWithConstruction(RBConstruction& rb_sys_in)
  : 
  rb_sys(rb_sys_in)
  {
  }
	
  /**
   * The RBConstruction object that will
   * use this assembly.
   */
  RBConstruction& rb_sys;
};

}

#endif // LIBMESH_ELEM_ASSEMBLY_H
