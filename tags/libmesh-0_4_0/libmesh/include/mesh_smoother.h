// $Id: mesh_smoother.h,v 1.2 2003-05-15 23:34:34 benkirk Exp $

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



#ifndef __mesh_smoother_h__
#define __mesh_smoother_h__



// C++ Includes   -----------------------------------

// forward declarations


// Local Includes -----------------------------------
#include "mesh.h"


/**
 * This class provides the necessary interface for
 * mesh smoothing.  Concrete mesh smoothers will be
 * derived from this abstract base class.
 *
 * \author John W. Peterson
 * \date 2002-2003
 * \version $Revision: 1.2 $
 */


// ------------------------------------------------------------
// MeshSmoother class definition
class MeshSmoother
{
public:
  /**
   * Constructor.  Sets the constant mesh reference
   * in the protected data section of the class.
   */
  MeshSmoother(Mesh& mesh) : _mesh(mesh) {}

  /**
   * Destructor.
   */
  virtual ~MeshSmoother() {}

  /**
   * Function which actually performs the smoothing operations.
   */
  virtual void smooth() = 0;
  
protected:

  Mesh& _mesh;
  
};


#endif
