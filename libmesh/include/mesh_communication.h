// $Id: mesh_communication.h,v 1.1 2003-02-28 23:37:43 benkirk Exp $

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



#ifndef __mesh_communication_h__
#define __mesh_communication_h__



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "mesh_common.h"


// Forward declarations
class Mesh;


/**
 * This is the \p MeshCommunication class.  It handles all the details
 * of communicating mesh information from one processor to another.  All
 * parallelization of the \p Mesh data structures is done via this class.
 *
 * @author Benjamin S. Kirk, 2003
 */


// ------------------------------------------------------------
// MeshCommunication class definition
class MeshCommunication
{
public:

  /**
   * Constructor.
   */
  MeshCommunication (Mesh& mesh);

  /**
   * Destructor.
   */
  ~MeshCommunication () {}

  /**
   * Clears all data structures and returns to a pristine state.
   */
  void clear ();
  
  /**
   * Finds all the processors that may contain
   * elements that neighbor my elements.  This list
   * is guaranteed to include all processors that border
   * any my elements, but may include additional ones as
   * well.
   */
  void find_neighboring_processors() {}

  
private:


  /**
   * Reference to the mesh.
   */
  Mesh& _mesh;
};



//--------------------------------------------------------------
// MeshCommunication inline members
inline
MeshCommunication::MeshCommunication (Mesh& mesh) :
  _mesh (mesh)
{
}



inline
void MeshCommunication::clear ()
{
}


#endif // end #ifndef __mesh_communication_h__ 
