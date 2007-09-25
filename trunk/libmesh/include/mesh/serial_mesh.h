// $Id: serial_mesh.h,v 1.1 2007-09-25 19:59:23 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __serial_mesh_h__
#define __serial_mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "unstructured_mesh.h"

/**
 * The \p SerialMesh class is derived from the \p MeshBase class,
 * and currently represents the default Mesh implementation.
 * Most methods for this class are found in MeshBase, and most
 * implementation details are found in UnstructuredMesh.
*/

// ------------------------------------------------------------
// Mesh class definition
class SerialMesh : public UnstructuredMesh
{
 public:

  /**
   * Constructor.  Requires the dimension and optionally
   * a processor id.  Note that \p proc_id should always
   * be provided for multiprocessor applications.
   */
  SerialMesh (unsigned int d);

  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  SerialMesh (const UnstructuredMesh& other_mesh) : 
    UnstructuredMesh(other_mesh) {}

  /**
   * Destructor.
   */
  ~SerialMesh();
};





#endif
