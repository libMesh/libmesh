// $Id: mesh_init.h,v 1.1 2003-02-10 03:55:51 benkirk Exp $

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



#ifndef __mesh_init_h__
#define __mesh_init_h__


// C++ includes

// Local includes
#include "mesh_common.h"


namespace libMesh
{

  /**
   * Initialize the library for use.  This will call
   * PetscInitialize if PETSC is available.  You must call
   * this method before using any of the library functionality.
   */
  void init (int argc=0, char** argv=NULL);

  /**
   * Stop using the mesh library.  This will call PetscFinalize()
   * if PETSC is available.
   */
  void close ();
};


#endif // #define __mesh_init_h__
